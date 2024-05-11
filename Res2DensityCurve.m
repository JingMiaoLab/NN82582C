function DensityCurve=Res2DensityCurve(res)
% This function accepts an AQUA project file 'res' and outputs a Inflares density curve as shown in Fig. 3c.
[xx,yy,~]=size(res.datOrg);
figure(1)
h_im = imagesc(mean(res.datOrg(:,:,1:20),3));colormap hot ; axis off
h=imellipse;
position=wait(h);
BW = createMask(h,h_im);
close(1)
[y,x]=ind2sub([xx,yy],find(BW==1));
y=round(mean(y));
x=round(mean(x));
DistMap=zeros(xx,yy);
for A=1:xx
    for B=1:yy
        DistMap(A,B)=dist([x,y],[A;B]);
    end
end
for A=1:length(res.ftsFav.loc.x2D)
    [X,Y]=ind2sub([xx,yy],res.ftsFav.loc.x2D{A});
    X=mean(X);Y=mean(Y);
    Dstdata(A)=dist([X,Y],[x;y]);
end
DensityCurve=zeros(1,300);
for A=1: 300
     % The default spacing between circular rings is 50 pixels.
    Region=(DistMap>=A&DistMap<A+50);
    DensityCurve(A)=length(find(Dstdata>=A&Dstdata<A+50))./(sum(Region(:)));
end
end