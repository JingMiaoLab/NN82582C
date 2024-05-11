function ROImap=Res2RoiRegion(res)
%This function takes an AQUA project file 'res' as input and 
%outputs a two-dimensional ROI of Inflares as shown in Fig. 2c, 
%with the ROI color mapping corresponding to the maximum response for the signal.
[x,y,~]=size(res.datOrg);
ROI2D=zeros(x,y);
for i=1:length(res.ftsFav.curve.tEnd)
ROI2D(res.ftsFav.loc.x2D{1,i})=res.ftsFav.curve.dffMax(1,i);
end
end


