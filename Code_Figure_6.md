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


### heatmap

```py
import glob
import pandas as pd
import numpy as np
import seaborn as sns

filelist=glob.glob("./result/03_diffexp/*res.tsv")
rawres=pd.DataFrame()
for i in filelist:
    tmp=pd.read_csv(i,sep='\t')
    rawres=rawres.append(tmp)

DEGres=rawres[(rawres.padj<0.05)&(rawres.log2FoldChange.abs()>1)]
DEgene=DEGres.drop_duplicates("gene")
DEgene.to_csv("DEA.genelist.csv",index=False)

info2=pd.read_csv("mm10.anno.tsv",sep='\t')
info2[info2.gene_name.isin(info1.gene_name)].to_csv("inflammatory_info.csv",index=False)

iDEG=DEgene.merge(info2[info2.gene_name.isin(d4.index)].drop("gene_name",axis=1),on='gene').drop(['chr','start','end','GOdes','GOtypes','GOnames','GOids','ENSEMBL'],axis=1)
iDEG=iDEG.loc[iDEG.log2FoldChange.abs().sort_values(ascending=False).index].reset_index(drop=True)
tmp=iDEG.applymap(lambda x:"%.3e"%x if isinstance(x,float) else x).drop("gene",axis=1)

tmp='''
ENSMUSG00000018774.13	Cd68
ENSMUSG00000040552.8	C3ar1
ENSMUSG00000021624.9	Cd180
ENSMUSG00000002603.15	Tgfb1
ENSMUSG00000030341.17	Tnfrsf1a
ENSMUSG00000025017.10	Pik3ap1
ENSMUSG00000027995.10	Tlr2
ENSMUSG00000030748.9	Il4ra
ENSMUSG00000019122.8	Ccl9
ENSMUSG00000052336.7	Cx3cr1
ENSMUSG00000018927.3	Ccl6
ENSMUSG00000020676.2	Ccl11
ENSMUSG00000035385.5	Ccl2
ENSMUSG00000028599.10	Tnfrsf1b
ENSMUSG00000024810.16	Il33
ENSMUSG00000021025.8	Nfkbia
ENSMUSG00000035373.2	Ccl7
ENSMUSG00000026072.12	Il1r1
ENSMUSG00000025888.6	Casp1
ENSMUSG00000008845.9	Cd163
ENSMUSG00000024164.15	C3
'''
from io import StringIO
tmp2=pd.read_csv(StringIO(tmp),sep='\t',header=None)
tmp3=iDEG[iDEG.gene.isin(tmp2.iloc[:,0])]

tmp=d4.loc[iDEG.gene_name]
labels=pd.Series(tmp.index.values).apply(lambda x:x if x in tmp3.gene_name.values else "").to_list()


sns.set(font_scale=0.8) 
sns.clustermap(tmp,cmap='coolwarm',cbar_pos=(1, .65, .03, .2),figsize=(10,10),yticklabels=labels,
dendrogram_ratio=[0.2,0.001],col_cluster=False,row_cluster=False,z_score=0)
```

### GSEA
```py
filelist=glob.glob("./result/03_diffexp/*res.tsv")
rawres=pd.DataFrame()
for i in filelist:
    tmp=pd.read_csv(i,sep='\t')
    rawres=rawres.append(tmp)

tpm=pd.read_csv("./result/02_quant/htseq.quant.tpm.txt",sep='\t',index_col=0)
tpm2=tpm[tpm.mean(axis=1)>1]
tpm2.shape

rawres2=rawres[rawres.gene.isin(tpm2.index.values)].dropna()
rawres2.contrast.value_counts()


gdata=rawres2[rawres2.contrast=="Condition_MCAO_vs_PB"].copy().dropna()
gdata.index=gdata.gene_name
geneList1=gdata.loc[:,'log2FoldChange'].sort_values(ascending=False)
geneList1.shape


gdata=rawres2[rawres2.contrast=="Condition_MCAO_vs_Sham"].copy().dropna()
gdata.index=gdata.gene_name
geneList2=gdata.loc[:,'log2FoldChange'].sort_values(ascending=False)
geneList2.shape

```


```r
geneList=sort(unlist(geneList1),decreasing = T)
orgDb=org.Mm.eg.db
gseares1 <- gseGO(geneList         = geneList,
                    OrgDb         = orgDb,
                    keyType       = 'SYMBOL',
                    ont='ALL')
gseadf1<-gseares1[order(gseares1$enrichmentScore, decreasing = T),]
saveRDS(gseares1,'Condition_MCAO_vs_PB.gsea.go.result0228.rds')


gseares1=readRDS("Condition_MCAO_vs_PB.gsea.go.result0228.rds")
pdf(file="Condition_MCAO_vs_PB.gsea.go.plot.pdf")
gseaplot2(gseares1,"GO:0006954",title="MCAO_vs_PB:inflammatory response")
dev.off()


geneList=sort(unlist(geneList2),decreasing = T)
orgDb=org.Mm.eg.db
gseares2 <- gseGO(geneList         = geneList,
                    OrgDb         = orgDb,
                    keyType       = 'SYMBOL',
                    ont='ALL',pvalueCutoff = 0.9)
gseadf2<-gseares2[order(gseares2$enrichmentScore, decreasing = T),]
saveRDS(gseares2,'Condition_MCAO_vs_Sham.gsea.go.result0228.rds')

pdf(file="Condition_MCAO_vs_Sham.gsea.go.plot.pdf")
gseares2=readRDS("Condition_MCAO_vs_Sham.gsea.go.result0228.rds")
gseaplot2(gseares2,"GO:0006954",title="MCAO_vs_Sham:inflammatory response")
dev.off()

```

