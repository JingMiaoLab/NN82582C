#  code.md

### run DESeq2

```r
library(tximport)
library(readr)
library(stringr)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(rtracklayer)


meta$Condition <- factor(meta$Condition)
meta$Group <- factor(meta$Group)

meta=meta%>%
    dplyr::filter(Condition%in%contrast)
rownames(meta) <- meta$Sample

samples<-meta%>%
    tibble::as_tibble()%>%
    dplyr::filter(Condition%in%contrast)%>%
    as.data.frame()

name_tmp=gsub("-",'.',samples$Sample)
rownames(samples)=gsub("\\+",'.',name_tmp)

df=df[,rownames(samples)]
df=df[rowSums(is.na(df))<(dim(df)[2]*0.8),]
cts=df[rowSums(df<2)<(dim(df)[2]*0.8),]
#cts[rowSums(cts)<10,]=0
#cts[rowSums(cts<2)>(dim(cts)[2]*0.8),]=0

if (length(unique(samples$Group))==1){
  ddSet <- DESeqDataSetFromMatrix(countData  = cts,
                                  colData = samples,
                                  design= ~ Condition)
}else{
  ddSet <- DESeqDataSetFromMatrix(countData  = cts,
                                  colData = samples,
                                  design= ~ Group + Condition)

}

ref=contrast[1]
ddSet$Condition <- relevel(ddSet$Condition, ref = ref)
if (length(contrast)>2){
    dds <- DESeq(ddSet, parallel=parallel)
}else{
    dds <- DESeq(ddSet,parallel=TRUE,test='LRT',reduced=~1)
}


geneMaptmp <- gtfinfo %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(type=="gene") 
if ('gene_name' %in% colnames(geneMaptmp)){
  geneMap<-geneMaptmp%>%
          dplyr::select(c(gene_name,gene_id)) 
}else if ('gene_symbol' %in% colnames(geneMaptmp)){
  geneMap<-geneMaptmp%>%
          dplyr::select(c(gene_symbol,gene_id))
  colnames(geneMap)<-c("gene_name",'gene_id')
} else{
  geneMap<-geneMaptmp%>%
        dplyr::select(c(gene_id))
  geneMap$gene_name=geneMap$gene_id
}


id_trans <-function(data_input){
  data_input <-data_input%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var='gene') %>% 
    dplyr::left_join(geneMap,by=c('gene'='gene_id'))%>% 
    dplyr::select(gene,gene_name, everything())
} 

get_result<-function(dds,resultname){
  #res <- results(dds, contrast=contrast,  parallel=parallel)
  res <- lfcShrink(dds, coef=resultname, type="apeglm",parallel=parallel)
  res%>%
    id_trans()%>%
    mutate(contrast=resultname)
}

if (length(unique(samples$Condition))>2){
    test<-'LRT'
}else{
    test<-'Wald'
}

Last=length(resultsNames(dds))
name=resultsNames(dds)[Last]
res=get_result(dds,name)

if (test=='LRT'){
  for (i in 3:length(resultsNames(dds))){
    name=resultsNames(dds)[i]
    tmp=get_result(dds,name)
    res<-res%>%
        bind_rows(tmp)
  }
}

if (dim(colData(dds))[1]<30){
    ndds <- rlog(dds, blind=FALSE)
} else {
    ndds <- vst(dds, blind=FALSE)
}
exp <- as.data.frame(assay(ndds)) %>% 
  id_trans() %>% 
  type_convert()

DESeq2_res_file=paste0(outdir,Name,'.DESeq2_res.tsv')
exp_file=paste0(outdir,Name,'.DESeq2_exp.tsv')
meta_file=paste0(outdir,Name,'.DESeq2_meta.tsv')

write_tsv(samples,meta_file)
write_tsv(exp,exp_file)
write_tsv(res,DESeq2_res_file)
```


### enrich

```py
import glob
import pandas as pd
import numpy as np
import seaborn as sns

result=pd.DataFrame()
for i in glob.glob("./result/03_diffexp/*DESeq2_res.tsv"):
    tmp=df=pd.read_csv(i,sep='\t')
    result=result.append(tmp)

gdata=result[(result.contrast=="Condition_Saline_FLA1_vs_PBN_FLA1")].copy().dropna()
gdata.index=gdata.gene_name
geneList1=gdata.loc[:,'log2FoldChange'].sort_values(ascending=False)


gdata=result[(result.contrast=="Condition_wt_vs_cko")].copy().dropna()
gdata.index=gdata.gene_name
geneList2=gdata.loc[:,'log2FoldChange'].sort_values(ascending=False)
```

```r

geneList=sort(unlist(geneList1),decreasing = T)
orgDb=org.Mm.eg.db
gseares1 <- gseGO(geneList         = geneList,
                    OrgDb         = orgDb,
                    keyType       = 'SYMBOL',
                    ont='ALL')

gseadf1<-gseares1[order(gseares1$enrichmentScore, decreasing = T),]
saveRDS(gseares1,'Condition_Saline_vs_PBN.gsea.go.result.rds')


geneList=sort(unlist(geneList2),decreasing = T)
orgDb=org.Mm.eg.db
gseares2 <- gseGO(geneList         = geneList,
                    OrgDb         = orgDb,
                    keyType       = 'SYMBOL',
                    ont='ALL')

gseadf2<-gseares2[order(gseares2$enrichmentScore, decreasing = T),]
saveRDS(gseares2,'Condition_wt_vs_ko.gsea.go.result.rds')


gseares1=readRDS("Condition_Saline_vs_PBN.gsea.go.result.rds")
gseaplot(gseares1,"GO:0050727",title="regulation of inflammatory response")



gseares2=readRDS("Condition_wt_vs_ko.gsea.go.result.rds")
gseaplot(gseares2,"GO:0050727",title="regulation of inflammatory response")
```


### gene plot
```py

def get_plot(DEname,grouplist,genelist=[]):
    data=res[res.contrast==DEname]
    if len(genelist)>0:
        data=data[data.gene_name.isin(genelist)]
    Meta=meta[meta.Group.isin(grouplist)]
    Quant=quant.loc[data.gene_name,Meta.Name].applymap(lambda x:np.log10(x+1))
    sdata=Quant.T.groupby(Meta.set_index("Name").Group).mean().T
    sdata=sdata.merge(data,left_index=True,right_on="gene_name")
    sdata.loc[:,"regulated"]=sdata.apply(lambda x:"0" if x.log2FoldChange>0 else "1",axis=1)
    sdata.loc[:,'FC']=sdata.log2FoldChange.abs()
    sdata=sdata.sort_values("FC",ascending=False)
    #genlist=sdata.groupby("regulated").head(10).gene_name
    genlist=sdata[sdata.regulated=='0'].groupby("regulated").head(10).gene_name
    genlist=genlist.append(pd.Series(['Il1b', 'Il6', 'Ccl2', 'Ccl5', 'Cxcl10', 'Tnf', 'Nfkb1']))
    sdata.loc[:,'Show']=sdata.gene_name.apply(lambda x :True if x in genlist.values else False)
    return sdata

info=pd.read_csv("mm10.anno.tsv",sep='\t')
genelist1=info[info.GOnames.fillna("").str.contains("inflammatory|proliferation|chemokine|cytokine|apoptosis")][['gene_name','GOnames']]

res=pd.DataFrame()
for i in glob.glob("./result/03_diffexp/*DEG_res.tsv"):
    tmp=df=pd.read_csv(i,sep='\t')
    res=res.append(tmp)

data2=get_plot("Condition_Saline_FLA1_vs_PBN_FLA1",['Saline','PBN'],genelist1.gene_name.to_list())  
data4=get_plot("Condition_wt_vs_cko",['wt','cko'],genelist1.gene_name.to_list())  
```


```r

p<-ggplot(data2, aes(PBN,Saline,color=regulated,label=ifelse(Show, gene_name, ""))) +
    geom_point()+
    scale_color_manual(values =c("red",'blue'))+
    geom_text_repel()+
    xlim(0,3.1)+ylim(0,3.1)+
    #scale_y_log10(limits = c(0.1, 1.5e3)) +scale_x_log10(limits = c(0.1, 1.5e3))+ 
    theme(legend.position = "none",panel.background = element_rect(fill = "white", colour = "grey50"),text = element_text(size = 20))+ 
    labs(title = "log10(TPM+1)")
    
ggsave(file="fig1.svg",plot=p)



p<-ggplot(data4, aes(cko,wt,color=regulated,label=ifelse(Show, gene_name, ""))) +
    geom_point()+
    scale_color_manual(values =c("red",'blue'))+
    geom_text_repel()+
    xlim(0,3.1)+ylim(0,3.1)+
    theme(legend.position = "none",panel.background = element_rect(fill = "white", colour = "grey50"),text = element_text(size = 20))+
    labs(title = "log10(TPM+1)")
ggsave("fig2.svg")
```