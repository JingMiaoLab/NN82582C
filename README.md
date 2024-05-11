# NN-A82582C
Main code for article NN-A82582C
The main script includes two parts: 

**1 Analysis script for image data**; 

**2 Analysis script for sequencing data**.

If you have any feedback or issue, you are welcome to either post issue in Issues section or send email to jingmiao@cibr.ac.cn.

# 图像数据
对于双光子成像数据[示例数据.tif](\ABCD.tif)的数据经过了 [AQuA](https://github.com/yu-lab-vt/AQuA)（xuelong mi et.al natureneuroscience 2019）进行了预处理,预设参数可以在 [此处](\预设.csv) 下载.
基于AQuA识别后的工程文件，请在[此处](\ABCD.mat)下载
1. 通过AQUA识别结果得到inflare ROI
2. 通过AQUA识别结果得到Inflare Density map
3. 通过AQUA识别结果得到Density Curve
4. 通过AQUA识别结果得到Inflare Frequency

# 测序数据
测序数据已在GEO公开，[Seq Data in Figure9](ww.1234), [Seq Data in Extend data Figure9](ww.1234)
