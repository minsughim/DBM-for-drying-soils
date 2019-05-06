day = 5;
%timeLine = timeLine(1:1440);

g=graph('2d-deg6-trilattice',{n,m,Lp});
listV = g.vertices;
listE = g.edges;
%periodic boundary condition for y direction
listE(m*n-1,m) = 1;
listE(m*n,m) = 1;
listE(m,m*n-1) = 1;
listE(m,m*n) = 1;
for i = 1:2:(m-1)
    listE((n-1)*m+i ,i) = 1;
    listE(i,(n-1)*m+i) = 1;
end
for i = 2:2:(m-2)
    listE((n-1)*m+i-1 ,i) = 1;
    listE((n-1)*m+i ,i) = 1;
    listE((n-1)*m+i+1 ,i) = 1;
    listE(i,(n-1)*m+i-1) = 1;
    listE(i,(n-1)*m+i) = 1;
    listE(i,(n-1)*m+i+1) = 1;
end

xlist = reshape(listV(:,1), [m,n])*1000;
ylist = reshape(listV(:,2), [m,n])*(-1000);


avgO2 = zeros(length(depthList),length(timeLine));
stdO2 = zeros(length(depthList),length(timeLine));
avgPH = zeros(length(depthList),length(timeLine));
stdPH = zeros(length(depthList),length(timeLine));
avgO2(:,:) = nanmean(timeConcDist{1},2);
avgPH(:,:) = nanmean(timeConcDist2{5},2);
stdO2(:,:) = nanstd(timeConcDist{1},0,2);
stdPH(:,:) = nanstd(timeConcDist2{5},0,2);

avgCO2 = zeros(length(depthList),length(timeLine));
stdCO2 = zeros(length(depthList),length(timeLine));
avgCO2(:,:) = nanmean(timeConcDist{2},2);
stdCO2(:,:) = nanstd(timeConcDist{2},0,2);

avgCH2O = zeros(length(depthList),length(timeLine));
stdCH2O = zeros(length(depthList),length(timeLine));
avgCH2O(:,:) = nanmean(timeConcDist{4},2);
stdCH2O(:,:) = nanstd(timeConcDist{4},0,2);

avgNH4 = zeros(length(depthList),length(timeLine));
stdNH4 = zeros(length(depthList),length(timeLine));
avgNH4(:,:) = nanmean(timeConcDist{6},2);
stdNH4(:,:) = nanstd(timeConcDist{6},0,2);


figure(1)
subplot(1,5,1)
mesh(timeLine,depthList, avgO2)
subplot(1,5,2)
mesh(timeLine,depthList, avgPH)
subplot(1,5,3)
mesh(timeLine,depthList, avgCO2)
subplot(1,5,4)
mesh(timeLine,depthList, avgCH2O)
subplot(1,5,5)
mesh(timeLine,depthList, avgNH4)

% 
% HoriList = 0:Lp:Lp*(n-1);
% 
% figure(1)
% timeT = 72*6-1; %before sunset
% temp = zeros(m,n);
% temp(:,:) = timeConcDist{1}(:,:,timeT);
% minDayO2avg = nanmean(temp,2);
% minDayO2std = nanstd(temp,0,2);
% 
% %temp(end,end) = max(timeConcDist{1}(:));
% temp(end,end) = 10;
% contour(HoriList*1000,depthList*1000,temp,0:0.1:10)
% figure(2)
% timeT = 72*8-1; %before sunrise
% temp = zeros(m,n);
% temp(:,:) = timeConcDist{1}(:,:,timeT);
% midNightO2avg = nanmean(temp,2);
% minNightO2std = nanstd(temp,0,2);
% %temp(end,end) = max(timeConcDist{1}(:));
% temp(end,end) = 10;
% contour(HoriList*1000,depthList*1000,temp,0:0.1:10)
% 
% 
% figure(3) 
% timeT = 72*6-1; %before sunset
% temp = zeros(m,n);
% temp(:,:) = timeConcDist2{5}(:,:,timeT);
% minDaypHavg = nanmean(temp,2);
% minDaypHstd = nanstd(temp,0,2);
% temp(end,end) = max(timeConcDist{1}(:));
% temp(end,end) = 11;
% contour(HoriList*1000,depthList*1000,temp,6:0.1:11)
% figure(4)
% timeT = 72*8-1; %before sunrise
% temp = zeros(m,n);
% temp(:,:) = timeConcDist2{5}(:,:,timeT);
% midNightpHavg = nanmean(temp,2);
% minNightpHstd = nanstd(temp,0,2);
% temp(end,end) = max(timeConcDist{1}(:));
% temp(end,end) = 11;
% contour(HoriList*1000,depthList*1000,temp,6:0.1:11)



Rad3Over2 = sqrt(3)/2;
Rad3 = sqrt(3);
LayerDepth = 0.02;
depthList = linspace(0,LayerDepth,m);
Lp = (depthList(2)-depthList(1))/Rad3Over2; % size of patch for toplayer [m]
dl = Lp/Rad3; % length of the edge of hexagon (top layer)
R1 = 10^(-4);
R2 = 10^(-7);
patchArea = Rad3Over2*dl*dl;

nightT = 1368;
dayT = 1224;

HoriList = 0:Lp:Lp*(n-1);
figure(2)
subplot(2,9,1)
contour(xlist,ylist, waterVolume/patchArea)
colorbar;
subplot(2,9,2)
contour(xlist,ylist, timeConcDist{1}(:,:,dayT))
colorbar;
subplot(2,9,3)
contour(xlist,ylist,timeConcDist{2}(:,:,dayT))
colorbar;
subplot(2,9,4)
contour(xlist,ylist, timeConcDist2{5}(:,:,dayT))
colorbar;
subplot(2,9,5)
contour(xlist,ylist, timeConcDist{6}(:,:,dayT)+timeConcDist2{3}(:,:,dayT))
colorbar;
subplot(2,9,6)
contour(xlist,ylist,timeConcDist{7}(:,:,dayT))
colorbar;
subplot(2,9,7)
contour(xlist,ylist, timeConcDist{5}(:,:,dayT))
colorbar;
subplot(2,9,8)
contour(xlist,ylist, timeConcDist{8}(:,:,dayT))
colorbar;
subplot(2,9,9)
contour(xlist,ylist, timeConcDist{4}(:,:,dayT))
colorbar;
subplot(2,9,10)
contour(xlist,ylist, waterVolume/patchArea)
colorbar;
subplot(2,9,11)
contour(xlist,ylist, timeConcDist{1}(:,:,nightT))
colorbar;
subplot(2,9,12)
contour(xlist,ylist,timeConcDist{2}(:,:,nightT))
colorbar;
subplot(2,9,13)
contour(xlist,ylist, timeConcDist2{5}(:,:,nightT))
colorbar;
subplot(2,9,14)
%contour(xlist,ylist, sitesC{6}+sitesC2{3})
contour(xlist,ylist, timeConcDist{6}(:,:,nightT)+timeConcDist2{3}(:,:,nightT))
colorbar;
subplot(2,9,15)
contour(xlist,ylist,timeConcDist{7}(:,:,nightT))
colorbar;
subplot(2,9,16)
contour(xlist,ylist, timeConcDist{5}(:,:,nightT))
colorbar;
subplot(2,9,17)
contour(xlist,ylist, timeConcDist{8}(:,:,nightT))
colorbar;
subplot(2,9,18)
contour(xlist,ylist, timeConcDist{4}(:,:,nightT))
colorbar;


figure(3)

timeT = 288*day-145; %before sunset
subplot(2,5,1)
contour(xlist,ylist, (popMovie{1,timeT}+popMovie{2,timeT}+popMovie{3,timeT}+popMovie{4,timeT})*10)
colorbar;
subplot(2,5,2)
contour(xlist,ylist, popMovie{5,timeT}*10)
colorbar;
subplot(2,5,3)
contour(xlist,ylist,popMovie{6,timeT}*10)
colorbar;
subplot(2,5,4)
contour(xlist,ylist, popMovie{7,timeT}*10)
colorbar;
subplot(2,5,5)
contour(xlist,ylist, popMovie{8,timeT}*10)
colorbar;
timeT = 288*day-1; %before sunrise
subplot(2,5,6)
contour(xlist,ylist, (popMovie{1,timeT}+popMovie{2,timeT}+popMovie{3,timeT}+popMovie{4,timeT})*10)
colorbar;
subplot(2,5,7)
contour(xlist,ylist, popMovie{5,timeT}*10)
colorbar;
subplot(2,5,8)
contour(xlist,ylist, popMovie{6,timeT}*10)
colorbar;
subplot(2,5,9)
contour(xlist,ylist,popMovie{7,timeT}*10)
colorbar;
subplot(2,5,10)
contour(xlist,ylist,popMovie{8,timeT}*10)
colorbar;

figure(10)
timeT = 288*day-145; %before sunset
temp = (popMovie{1,timeT}+popMovie{2,timeT}+popMovie{3,timeT}+popMovie{4,timeT})*10;
avgPop = mean(temp,2);
minPop = min(temp');
maxPop = max(temp');
areaT(:,1) = minPop;
areaT(:,2) = (maxPop-minPop)';
area(ylist(:,1),areaT)
hold on
plot(ylist(:,1),avgPop)

for i = 5:8
temp = popMovie{i,timeT}*10;
avgPop = mean(temp,2);
minPop = min(temp');
maxPop = max(temp');
areaT(:,1) = minPop;
areaT(:,2) = (maxPop-minPop)';
area(ylist(:,1),areaT)
hold on
plot(ylist(:,1),avgPop)
end


figure(4)
%semilogy(timeLine/24, effluxList)
totVerticalArea = mean(totaldM(:))*Lp*n;
%semilogy(timeLine/24, -1*effluxList*3600/dt/plottt/totVerticalArea) %g/m^2hr
%hold on
%semilogy(timeLine/24, effluxList*3600/dt/plottt/totVerticalArea) %g/m^2hr
plot(timeLine/24, effluxList*10^6/dt/plottt/totVerticalArea) %\mug/m^2s
hold on
plot(timeLine/24, effluxList2*10^6/dt/plottt/totVerticalArea) %\mug/m^2s
%hold on
%plot(timeLine/24, effluxList3*10^6/dt/plottt/totVerticalArea) %\mug/m^2s

figure(5)
examineT1 = length(Mumean);
semilogy(timeLine(1:2:examineT1)/24, popWalkers(1:2:examineT1,:).*Mumean(1:2:examineT1,:))
hold on
semilogy(timeLine(2:2:examineT1)/24, popWalkers(2:2:examineT1,:).*Mumean(2:2:examineT1,:))

figure(6)
semilogy(timeLine(1:examineT1)/24, popWalkers(1:examineT1,:))
hold on
semilogy(timeLine(1:examineT1)/24, dormpopWalkers(1:examineT1,:),'--')

figure(7)
timeT = 288*day-145;
hold on
errorbar(depthList(1:50), avgPH((1:50),timeT), stdPH((1:50),timeT))
timeT = 288*day-1;
hold on
errorbar(depthList(1:50), avgPH((1:50),timeT), stdPH((1:50),timeT))


figure(8)
timeT = 288*day-145;
hold on
plot(avgPH(:,timeT),-1000*depthList)
timeT = 288*day-1;
hold on
plot(avgPH(:,timeT),-1000*depthList)

figure(9)
HoriList = 0:Lp:Lp*(n-1);
timeT = 288*day-145;
subplot(1,2,1)
hold on
temp = zeros(m,n);
temp(:,:) = timeConcDist2{5}(:,:,timeT);
temp(end,1) = 3;
temp(end,end) = 9;
contour(xlist,ylist,temp,100)
colorbar;
subplot(1,2,2)
hold on
timeT = 288*day-1;
temp = zeros(m,n);
temp(:,:) = timeConcDist2{5}(:,:,timeT);
temp(end,1) = 3;
temp(end,end) = 9;
contour(xlist,ylist,temp,100)
colorbar;

