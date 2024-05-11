# Spatiotemporally selective astrocytic ATP dynamics encode injury information sensed by microglia following brain injury in mice
The main code for article NN-A82582C
which includes two parts: 

**1 Analysis code for image data**; 

**2 Analysis code for sequencing data**.

If you have any feedback or issue, you are welcome to either post issue in Issues section or send email to jingmiao@cibr.ac.cn.

## ATP1.0 *in vivo* imaging data analysis

### Exampluar data

You can download the exampluar data [here](\ABCD.tif) 

You can download the project files post signal extraction via [AQuA](https://github.com/yu-lab-vt/AQuA) (*Yizhi Wang et.al nat. neurosci 2019*)  from [here](\ABCD.mat).

Opts for signal extraction can be downloaded [here](/OptsOfAqua.csv).

### Code 

1. [**Res2DensityMap**](/Res2DensityMap.m) This function accepts an AQUA project file 'res' and outputs a normalized Inflares density map as shown in Fig. 3a.
2. [**Res2DensityCurve**](/Res2DensityCurve.m) This function accepts an AQUA project file 'res' and outputs a Inflares density curve as shown in Fig. 3c.
3. [**Res2RoiRegion**](/Res2RoiRegion.m) This function takes an AQUA project file 'res' as input and outputs a two-dimensional ROI of Inflares as shown in Fig. 2c, with the ROI color mapping corresponding to the maximum response for the signal.
4. [**Res2FreqPlt**](/Res2FreqPlt.m) This function accepts an AQUA project file 'res' and plot frequency figure as shown in Fig. 2d.

## Sequencing data analysis
### Raw data
The sequencing data has been uploaded to NCBI and can be downloaded here. 

[Seq Data in Figure 6](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP501521)

[Seq Data in Extend data Figure 9](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP501522)


