examineDay =5;
examineDay2 = 3;
pCO2frac = 1;
NfixationR = 0.1;
amountCations = 0.1;
avgT = 25;
peneD = 0.2;
reducedLevelC = 0.1;
pot1 = 3;
plottt = 5;
dataPoints = 1440;
depthList = 100;
numberOfSamples = 1;
avgO2 = zeros(depthList,dataPoints,numberOfSamples);
stdO2 = zeros(depthList,dataPoints,numberOfSamples);
avgPH = zeros(depthList,dataPoints,numberOfSamples);
stdPH = zeros(depthList,dataPoints,numberOfSamples);
avgCH2O = zeros(depthList,dataPoints,numberOfSamples);
stdCH2O = zeros(depthList,dataPoints,numberOfSamples);
effluxListAll = zeros(dataPoints,5, numberOfSamples);
effluxListAll2 = zeros(dataPoints,5, numberOfSamples);
timeLinelist = zeros(dataPoints,5);
MumeanAll = zeros(floor(dataPoints*0.5),8, numberOfSamples);


for indexS = 1:numberOfSamples

    serial_id = sprintf('HT_BSC_noPG_nodeposit_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);

    %serial_id = sprintf('HT_BSC_noPG_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);
    %serial_id = sprintf('HT_BSC_Day%dPot%.1fredC%.2findex%d',examineDay,pot1,reducedLevelC, indexS);
    %serial_id = sprintf('HT_BSC_intraPG2_Day%dPot%.1fredC%.2findex%d',examineDay,pot1,reducedLevelC, indexS);
    cd(serial_id);
   % serial_idtemp = sprintf('Dinural_cycle');
   % cd(serial_idtemp);
    load('Final.mat', 'listV', 'm' ,'n','timeConcDist','timeConcDist2','totaldM','Lp','porosity','dt','plottt','effluxList','popWalkers','Mumean','examineT','timeLine', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesC','sitesC2','sitesCg')
   % cd ..
    cd ..
    xlist = reshape(listV(:,1), [m,n])*1000;
    ylist = reshape(listV(:,2), [m,n])*(-1000);
    avgO2(:,:,indexS) = nanmean(timeConcDist{1},2);
    avgPH(:,:,indexS) = nanmean(timeConcDist2{5},2);
    stdO2(:,:,indexS) = nanstd(timeConcDist{1},0,2);
    stdPH(:,:,indexS) = nanstd(timeConcDist2{5},0,2);
    avgCH2O(:,:,indexS) = nanmean(timeConcDist{4},2);
    stdCH2O(:,:,indexS) = nanstd(timeConcDist{4},0,2);
    
    
    [effluxList2, timeSeriesDeltaM] = EffluxCalculation(depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg);

    totVerticalArea = mean(totaldM(:))*n*Lp;
    %totVerticalArea = sum(sum(totaldM(:,:).*InvasedIsland))*Lp;
   % effluxListAll2(:,:,indexS) = effluxList2*10^6/44/dt/plottt/totVerticalArea; %mumol/m^s sec
    effluxListAll(:,:,indexS) = effluxList*10^6/44/dt/plottt/totVerticalArea; %mumol/m^s sec
    
    examineT1 = length(timeLine);
    MumeanAll(:,:, indexS) = popWalkers(1:2:examineT1,:).*Mumean(1:2:examineT1,:);
    timeLinelist(:,indexS) = timeLine;
end

figure(1)
hold on
temp = mean(effluxListAll,3);
temp2 = std(effluxListAll,0,3);
errorbar(timeLine(1:dataPoints)/24, temp(:,2), temp2(:,2))
hold on
plot(timeLine(1:dataPoints)/24, temp(:,2))
hold on
%temp3 = mean(effluxListAll2,3);
%temp4 = std(effluxListAll2,0,3);
%errorbar(timeLinelist(:,2)/24, temp3(:,2), temp4(:,2))
%hold on
%plot(timeLine/24, temp3(:,2))
hold on
load('Rajeev2013_data.mat')
plot(RajeevCO2(:,1)/24, RajeevCO2(:,7:11)*10^6/44/3600)
saveas(h,CO2efflux,'jpg')


figure(2)
mean(avgO2,3);
imagesc(timeLine/24,depthList,ans);
figure(3)
mean(avgPH,3);
imagesc(timeLine/24,depthList,ans);
figure(4)
mean(avgCH2O,3);
imagesc(timeLine/24,depthList,ans);


meanO2profile = mean(avgO2,3);
stdO2profile = std(avgO2,0,3);
meanPHprofile = mean(avgPH,3);
stdPHprofile = std(avgPH,0,3);



O2plots = zeros(100,5);
O2plots(:,1) = depthList;
O2plots(:,2) = meanO2profile(:,1224);
O2plots(:,3) = stdO2profile(:,1224);
O2plots(:,4) = meanO2profile(:,1368);
O2plots(:,5) = stdO2profile(:,1368);
%save(strcat('./../../Fig/O2plots.txt'),'O2plots','-ASCII')

meanPHprofile = mean(avgPH,3);
stdPHprofile = std(avgPH,0,3);
PHplots = zeros(100,5);
PHplots(:,1) = depthList;
PHplots(:,2) = meanPHprofile(:,1224);
PHplots(:,3) = stdPHprofile(:,1224);
PHplots(:,4) = meanPHprofile(:,1368);
PHplots(:,5) = stdPHprofile(:,1368);
%save(strcat('./../../Fig/PHplots.txt'),'PHplots','-ASCII')
