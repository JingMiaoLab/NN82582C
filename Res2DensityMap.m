function dst = Res2DensityMap(res, fw)
% This function receives the engineering file 'res' identified by AQUA; 
% based on the Signal coordinate information provided therein, 
% it uses the Kernel Density Estimation Toolbox for MATLAB to calculate the norm inflare density map as shown in Figure 3a.
% Kernel Density Estimation Toolbox: https://ics.uci.edu/~ihler/code/kde.html
if nargin < 2
    fw = 50; %fw is used to control the range of the loss function, with a default value of 50
end
[xx, yy, ~] = size(res.datOrg);
for A = 1:length(res.ftsFav.curve.tEnd)
    [y, x] = ind2sub([xx, yy], res.ftsFav.loc.x2D{A});
    x = mean(x); y = mean(y);
    po(A, :) = [x, y];
end
po = [po; [0, 0; xx, yy]];
cd('E:\') %Need to change the Kde path
p = kde(po', [fw; fw], [], 'G');
dst = hist(p, [xx, yy]);
dst=dst./max(dst(:));
end
