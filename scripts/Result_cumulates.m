numberOfSamples = 5;
VariableShift = zeros(numberOfSamples,9);

VariableShift(1,:) = [5 3 1 0.1 0.1 25 0.2 0.1 1];
VariableShift(2,:) = [5 3 1 0.1 0.1 25 0.2 0.1 2];
VariableShift(3,:) = [5 3 1 0.1 0.5 25 0.2 0.1 1];
VariableShift(4,:) = [5 3 1 0.1 0.5 25 0.2 0.1 2];
VariableShift(5,:) = [5 3 2 0.1 0.1 25 0.2 0.1 1];
%VariableShift(6,:) = [5 3 2 0.1 0.1 25 0.2 0.1 2];


for indexT = 1:numberOfSamples
    
    examineDay = VariableShift(indexT,1);
    pot1 = VariableShift(indexT,2);
    pCO2frac = VariableShift(indexT,3);
    NfixationR = VariableShift(indexT,4);
    amountCations = VariableShift(indexT,5);
    avgT = VariableShift(indexT,6);
    peneD =VariableShift(indexT,7);
    reducedLevelC = VariableShift(indexT,8);
    indexS = VariableShift(indexT,9);
    
    %serial_id = sprintf('HT_BSC_noPG_nodeposit_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);
    serial_id = sprintf('HT_BSC_noPG_nodeposit_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);
    cd(serial_id);
    % serial_idtemp = sprintf('Dinural_cycle');
    % cd(serial_idtemp);
    load('Final.mat', 'listV', 'm' ,'n','timeConcDist','timeConcDist2','totaldM','Lp','porosity','dt','plottt','effluxList','popWalkers','Mumean','examineT','timeLine', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesC','sitesC2','sitesCg')
    % cd ..
    cd ..
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
    
    
    meanO2profile = mean(avgO2,3);
    stdO2profile = std(avgO2,0,3);
    O2plots = zeros(100,5);
    O2plots(:,1) = depthList;
    O2plots(:,2) = meanO2profile(:,1224);
    O2plots(:,3) = stdO2profile(:,1224);
    O2plots(:,4) = meanO2profile(:,1368);
    O2plots(:,5) = stdO2profile(:,1368);
    
    meanPHprofile = mean(avgPH,3);
    stdPHprofile = std(avgPH,0,3);
    PHplots = zeros(100,5);
    PHplots(:,1) = depthList;
    PHplots(:,2) = meanPHprofile(:,1224);
    PHplots(:,3) = stdPHprofile(:,1224);
    PHplots(:,4) = meanPHprofile(:,1368);
    PHplots(:,5) = stdPHprofile(:,1368);
    
    save(strcat('./',serial_id,'/','O2plots.txt'),'O2plots','-ASCII'); %plots midday and midnight
    save(strcat('./',serial_id,'/','PHplots.txt'),'PHplots','-ASCII');
    
    O2DistDay = zeros(m,n);
    O2DistNight = zeros(m,n);
    O2DistDay = timeConcDist{1}(:,:,1224);
    O2DistNight = timeConcDist{1}(:,:,1368);
    
    save(strcat('./',serial_id,'/','O2DistDay.txt'),'O2DistDay','-ASCII');
    save(strcat('./',serial_id,'/','O2DistNight.txt'),'O2DistNight','-ASCII');
    
    PHDistDay = zeros(m,n);
    PHDistNight = zeros(m,n);
    PHDistDay = timeConcDist{1}(:,:,1224);
    PHDistNight = timeConcDist{1}(:,:,1368);
    
    save(strcat('./',serial_id,'/','PHDistDay.txt'),'PHDistDay','-ASCII');
    save(strcat('./',serial_id,'/','PHDistNight.txt'),'PHDistNight','-ASCII');
    
    
    totVerticalArea = mean(totaldM(:))*n*Lp;
    effluxListAll = effluxList*10^6/44/dt/plottt/totVerticalArea; %mumol/m^s sec
    examineT1 = length(timeLine);
    MumeanAll = popWalkers(1:2:examineT1,:).*Mumean(1:2:examineT1,:);
    timeLinelist = timeLine;
    
    
    save(strcat('./',serial_id,'/','timeLine.txt'),'timeLine','-ASCII');
    save(strcat('./',serial_id,'/','MumeanAll.txt'),'MumeanAll','-ASCII');
    save(strcat('./',serial_id,'/','effluxListAll.txt'),'effluxListAll','-ASCII');
    save(strcat('./Results/',serial_id,'effluxList.txt'),'effluxListAll','-ASCII');
    save(strcat('./Results/',serial_id,'O2plots.txt'),'O2plots','-ASCII');
    save(strcat('./Results/',serial_id,'PHplots.txt'),'PHplots','-ASCII');    
end


VariableShift(1,:) = [5 3 1 0.1 0.1 25 0.2 0.1 1];
VariableShift(2,:) = [5 3 1 0.1 0.1 25 0.2 0.1 2];
VariableShift(3,:) = [5 3 1 0.1 0.5 25 0.2 0.1 1];
%VariableShift(4,:) = [5 3 1 0.1 0.5 25 0.2 0.1 2];
VariableShift(4,:) = [5 3 2 0.1 0.1 25 0.2 0.1 1];
%VariableShift(6,:) = [5 3 2 0.1 0.1 25 0.2 0.1 2];

for indexT = 1:numberOfSamples
    
    examineDay = VariableShift(indexT,1);
    pot1 = VariableShift(indexT,2);
    pCO2frac = VariableShift(indexT,3);
    NfixationR = VariableShift(indexT,4);
    amountCations = VariableShift(indexT,5);
    avgT = VariableShift(indexT,6);
    peneD =VariableShift(indexT,7);
    reducedLevelC = VariableShift(indexT,8);
    indexS = VariableShift(indexT,9);
    
    %serial_id = sprintf('HT_BSC_noPG_nodeposit_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);
    serial_id = sprintf('HT_BSC_noPG_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);
    cd(serial_id);
    % serial_idtemp = sprintf('Dinural_cycle');
    % cd(serial_idtemp);
    load('Final.mat', 'listV', 'm' ,'n','timeConcDist','timeConcDist2','totaldM','Lp','porosity','dt','plottt','effluxList','popWalkers','Mumean','examineT','timeLine', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesC','sitesC2','sitesCg')
    % cd ..
    cd ..
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
    
    
    meanO2profile = mean(avgO2,3);
    stdO2profile = std(avgO2,0,3);
    O2plots = zeros(100,5);
    O2plots(:,1) = depthList;
    O2plots(:,2) = meanO2profile(:,1224);
    O2plots(:,3) = stdO2profile(:,1224);
    O2plots(:,4) = meanO2profile(:,1368);
    O2plots(:,5) = stdO2profile(:,1368);
    
    meanPHprofile = mean(avgPH,3);
    stdPHprofile = std(avgPH,0,3);
    PHplots = zeros(100,5);
    PHplots(:,1) = depthList;
    PHplots(:,2) = meanPHprofile(:,1224);
    PHplots(:,3) = stdPHprofile(:,1224);
    PHplots(:,4) = meanPHprofile(:,1368);
    PHplots(:,5) = stdPHprofile(:,1368);
    
    save(strcat('./',serial_id,'/','O2plots.txt'),'O2plots','-ASCII'); %plots midday and midnight
    save(strcat('./',serial_id,'/','PHplots.txt'),'PHplots','-ASCII');
    
    O2DistDay = zeros(m,n);
    O2DistNight = zeros(m,n);
    O2DistDay = timeConcDist{1}(:,:,1224);
    O2DistNight = timeConcDist{1}(:,:,1368);
    
    save(strcat('./',serial_id,'/','O2DistDay.txt'),'O2DistDay','-ASCII');
    save(strcat('./',serial_id,'/','O2DistNight.txt'),'O2DistNight','-ASCII');
    
    PHDistDay = zeros(m,n);
    PHDistNight = zeros(m,n);
    PHDistDay = timeConcDist{1}(:,:,1224);
    PHDistNight = timeConcDist{1}(:,:,1368);
    
    save(strcat('./',serial_id,'/','PHDistDay.txt'),'PHDistDay','-ASCII');
    save(strcat('./',serial_id,'/','PHDistNight.txt'),'PHDistNight','-ASCII');
    
    
    totVerticalArea = mean(totaldM(:))*n*Lp;
    effluxListAll = effluxList*10^6/44/dt/plottt/totVerticalArea; %mumol/m^s sec
    examineT1 = length(timeLine);
    MumeanAll = popWalkers(1:2:examineT1,:).*Mumean(1:2:examineT1,:);
    timeLinelist = timeLine;
    
    
    save(strcat('./',serial_id,'/','timeLine.txt'),'timeLine','-ASCII');
    save(strcat('./',serial_id,'/','MumeanAll.txt'),'MumeanAll','-ASCII');
    save(strcat('./',serial_id,'/','effluxListAll.txt'),'effluxListAll','-ASCII');
    save(strcat('./Results/',serial_id,'effluxList.txt'),'effluxListAll','-ASCII');
    save(strcat('./Results/',serial_id,'O2plots.txt'),'O2plots','-ASCII');
    save(strcat('./Results/',serial_id,'PHplots.txt'),'PHplots','-ASCII');
    
end

%
% figure(1)
% hold on
% temp = mean(effluxListAll,3);
% temp2 = std(effluxListAll,0,3);
% errorbar(timeLine(1:dataPoints)/24, temp(:,2), temp2(:,2))
% hold on
% plot(timeLine(1:dataPoints)/24, temp(:,2))
% hold on
% %temp3 = mean(effluxListAll2,3);
% %temp4 = std(effluxListAll2,0,3);
% %errorbar(timeLinelist(:,2)/24, temp3(:,2), temp4(:,2))
% %hold on
% %plot(timeLine/24, temp3(:,2))
% hold on
% load('Rajeev2013_data.mat')
% plot(RajeevCO2(:,1)/24, RajeevCO2(:,7:11)*10^6/44/3600)
% saveas(h,CO2efflux,'jpg')
%
%
% figure(2)
% mean(avgO2,3);
% imagesc(timeLine/24,depthList,ans);
% figure(3)
% mean(avgPH,3);
% imagesc(timeLine/24,depthList,ans);
% figure(4)
% mean(avgCH2O,3);
% imagesc(timeLine/24,depthList,ans);
%
%
%
% meanO2profile = mean(avgO2,3);
% stdO2profile = std(avgO2,0,3);
% meanPHprofile = mean(avgPH,3);
% stdPHprofile = std(avgPH,0,3);
%
%
% O2plots = zeros(100,5);
% O2plots(:,1) = depthList;
% O2plots(:,2) = meanO2profile(:,1224);
% O2plots(:,3) = stdO2profile(:,1224);
% O2plots(:,4) = meanO2profile(:,1368);
% O2plots(:,5) = stdO2profile(:,1368);
% %save(strcat('./../../Fig/O2plots.txt'),'O2plots','-ASCII')
%
% meanPHprofile = mean(avgPH,3);
% stdPHprofile = std(avgPH,0,3);
% PHplots = zeros(100,5);
% PHplots(:,1) = depthList;
% PHplots(:,2) = meanPHprofile(:,1224);
% PHplots(:,3) = stdPHprofile(:,1224);
% PHplots(:,4) = meanPHprofile(:,1368);
% PHplots(:,5) = stdPHprofile(:,1368);


