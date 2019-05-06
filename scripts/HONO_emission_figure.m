 
load('/Users/mkim/Downloads/Constant_dark_water_lineardecay_Drying1_varyingT 1/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')
%load('/Users/mkim/Downloads/Constant_T_lineardecay_Drying1_varyingT/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')

figure(1)
hold on
plot(timeLine/24,totalNutri)

figure(2)
subplot(2,1,1)
totVerticalArea = mean(totaldM(:))*Lp*n;
hold on
plot(timeLine(1153:end)-4*24, effluxList2(1153:end,3)*14/17*10^9/dt/plottt/totVerticalArea) %\mug/m^2s
plot(timeLine(1153:end)-4*24, effluxList2(1153:end,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s
subplot(2,1,2)
hold on
plot(timeLine(1153:end)-4*24,avgWC(1153:end))

figure(3)
hold on
plot(avgWC, effluxList2*10^6/dt/plottt/totVerticalArea) %\mug/m^2s

figure(4)
hold on
plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s


load('/Users/mkim/Downloads/Constant_dark_water_lineardecay_Drying2_varyingT 1/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')
%load('/Users/mkim/Downloads/Constant_T_lineardecay_Drying2_varyingT/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')

figure(1)
hold on
plot(timeLine/24,totalNutri)

figure(2)
subplot(2,1,1)
totVerticalArea = mean(totaldM(:))*Lp*n;
hold on
plot(timeLine(1153:end)-4*24, effluxList2(1153:end,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

plot(timeLine(1153:end)-4*24, effluxList2(1153:end,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s
subplot(2,1,2)
hold on
plot(timeLine(1153:end)-4*24,avgWC(1153:end))
figure(3)
hold on
plot(avgWC, effluxList2*10^6/dt/plottt/totVerticalArea) %\mug/m^2s

figure(4)
hold on
plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s



load('/Users/mkim/Downloads/Constant_dark_water_lineardecay_Drying3_varyingT 1/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')
%load('/Users/mkim/Downloads/Constant_T_lineardecay_Drying3_varyingT/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')

figure(1)
hold on
plot(timeLine/24,totalNutri)

figure(2)
subplot(2,1,1)
totVerticalArea = mean(totaldM(:))*Lp*n;
hold on
plot(timeLine(1153:end)-4*24, effluxList2(1153:end,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

plot(timeLine(1153:end)-4*24, effluxList2(1153:end,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s
subplot(2,1,2)
hold on
plot(timeLine(1153:end)-4*24,avgWC(1153:end))
figure(3)
hold on
plot(avgWC, effluxList2*10^6/dt/plottt/totVerticalArea) %\mug/m^2s

figure(4)
hold on
plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s



load('/Users/mkim/Downloads/Constant_dark_water_lineardecay_Drying4_varyingT 1/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')
%load('/Users/mkim/Downloads/Constant_T_lineardecay_Drying4_varyingT/Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC')

figure(1)
hold on
plot(timeLine/24,totalNutri)
figure(2)
subplot(2,1,1)
totVerticalArea = mean(totaldM(:))*Lp*n;
hold on
plot(timeLine(1153:end)-4*24, effluxList2(1153:end,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

plot(timeLine(1153:end)-4*24, effluxList2(1153:end,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s
subplot(2,1,2)
hold on
plot(timeLine(1153:end)-4*24,avgWC(1153:end))
figure(3)
hold on
plot(avgWC, effluxList2*10^6/dt/plottt/totVerticalArea) %\mug/m^2s

figure(4)
hold on
plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,3)*17/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s
plot(avgWC(1153:1440)*100/0.24, effluxList2(1153:1440,4)*14/47*10^9/dt/plottt/totVerticalArea) %\mug/m^2s

