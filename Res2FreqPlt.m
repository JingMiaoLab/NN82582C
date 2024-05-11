function Res2FreqPlt(res)
figure
[~,~,tt]=size(res.datOrg);
for A=1:length(res.ftsFav.curve.tEnd)
    hold on
    X=[res.ftsFav.curve.tBegin(A),res.ftsFav.curve.tBegin(A),res.ftsFav.curve.tEnd(A),res.ftsFav.curve.tEnd(A)];
    Y=[-res.ftsFav.curve.dffMax(A),res.ftsFav.curve.dffMax(A),res.ftsFav.curve.dffMax(A),-res.ftsFav.curve.dffMax(A)];
    fill(X,Y,[0,0,0],'LineStyle','none')
    xlim([0,tt]) ;ylim([-2,2])
    hold off
end