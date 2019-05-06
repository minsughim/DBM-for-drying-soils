numberOfSamples = 8;
VariableShift = zeros(numberOfSamples,8);

desiccationIndex = 7;
hono = 0.005;
nh3 = 20;
%[0.1 5 10 15]
VariableShift(1,:) = [5 3 1 0.1 25 nh3 hono 1];
VariableShift(2,:) = [5 3 1 0.1 25 nh3 hono 2];
VariableShift(3,:) = [5 3 1 0.1 25 nh3 hono 3];
VariableShift(4,:) = [5 3 1 0.1 25 nh3 hono 4];
VariableShift(5,:) = [5 3 1 0.1 25 nh3 hono 5];
VariableShift(6,:) = [5 3 1 0.1 25 nh3 hono 6];
VariableShift(7,:) = [5 3 1 0.1 25 nh3 hono 7];
VariableShift(8,:) = [5 3 1 0.1 25 nh3 hono 8];
%VariableShift(8,:) = [5 3 1 0.1 25 3 8];
EffluxAll = zeros(1728,5,numberOfSamples);
totalNutriAll = zeros(1728,8,numberOfSamples);
avgWCAll = zeros(1728,numberOfSamples);
avgWFAll = zeros(1728,numberOfSamples);
alphaCO = 3;


waterContentsDynamicsTot = cell(numberOfSamples);
waterVolumeDynamicsTot = cell(numberOfSamples);
waterFilmDynamicsTot = cell(numberOfSamples);
InvasedIslandDynamicsTot = cell(numberOfSamples);
gasVolumeDynamicsTot = cell(numberOfSamples);


for indexT = 1:numberOfSamples
    examineDay = VariableShift(indexT,1);
    pot1 = VariableShift(indexT,2);
    pCO2frac = VariableShift(indexT,3);
    amountCations = VariableShift(indexT,4);
    avgT = VariableShift(indexT,5);
    NH3ppb = VariableShift(indexT,6);
    HONOppb =VariableShift(indexT,7);
    indexS = VariableShift(indexT,8);
    
    serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);
    cd(serial_id_new);
    serial_id3 = sprintf('WaterDynamics_Drying%d_varyingT',desiccationIndex);
    load(strcat(serial_id3,'.mat'));
    waterContentsDynamicsTot{indexT} = waterContentsDynamics;
    waterVolumeDynamicsTot{indexT} = waterVolumeDynamics;
    waterFilmDynamicsTot{indexT} = waterFilmDynamics;
    InvasedIslandDynamicsTot{indexT} = InvasedIslandDynamics;
    gasVolumeDynamicsTot{indexT} = gasVolumeDynamics;  
    serial_idtemp = sprintf('Constant_dark_nobio_H2ONO_NH3%d_HONO%d_Drying%d_varyingT', NH3ppb,HONOppb,desiccationIndex);
    cd(serial_idtemp);
    load('Final.mat','timeLine', 'timeConcDist', 'timeConcDist2', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC','avgWF','patchArea','TotalporeD','porosityM')
    cd ..
    cd ..
    totSoilVolume = patchArea*TotalporeD./porosityM;
    totSoilWeight = totSoilVolume*1300000;

    totalpH = zeros(1728,1);
    minpH = zeros(1728,1);
    maxpH = zeros(1728,1);
    meanpH = zeros(1728,1);
    
    for i = 1:length(timeLine)
       
        temp = waterVolumeDynamics(:,:,i).*timeConcDist2{4}(:,:,i);
        temp2 = waterVolumeDynamics(:,:,i);
        totalpH(i) = sum(temp(:))./sum(temp2(:));
        minpH(i) = min(min(timeConcDist2{5}(:,:,i)));
        maxpH(i) = max(max(timeConcDist2{5}(:,:,i)));
        meanpH(i) = mean(mean(timeConcDist2{5}(:,:,i)));
        
    end
    
    totVerticalArea = mean(totaldM(:))*Lp*n;
    EffluxAll(:,:,indexT) = effluxList2(:,:)*10^9/dt/plottt./totVerticalArea;
    totalNutriAll(:,:,indexT) = totalNutri(:,:);
    avgWCAll(:,indexT) = avgWC;
    avgWFAll(:,indexT) = avgWF;

    totalpHtot(:,indexT) = totalpH;
    minpHtot(:,indexT) = minpH;
    maxpHtot(:,indexT) = maxpH;
    meanpHtot(:,indexT) = meanpH;
       
end

serial_id2 = sprintf('HT_BSC_photoY3_Constant_dark_nobio_H2ONO_NH3%d_HONO%d_Drying%d',NH3ppb, HONOppb, desiccationIndex);
save(strcat('./Results_HONO_NH3/',serial_id2,'.mat'),'-v7.3');

