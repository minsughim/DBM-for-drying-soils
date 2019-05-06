numberOfSamples = 8;
VariableShift = zeros(numberOfSamples,9);

desiccationIndex = 4;
newT = 25;
hono = 1;
nh3 = 20;
%[0.1 5 10 15]
VariableShift(1,:) = [5 3 1 0.1 25 nh3 hono 1 newT];
VariableShift(2,:) = [5 3 1 0.1 25 nh3 hono 2 newT];
VariableShift(3,:) = [5 3 1 0.1 25 nh3 hono 3 newT];
VariableShift(4,:) = [5 3 1 0.1 25 nh3 hono 4 newT];
VariableShift(5,:) = [5 3 1 0.1 25 nh3 hono 5 newT];
VariableShift(6,:) = [5 3 1 0.1 25 nh3 hono 6 newT];
VariableShift(7,:) = [5 3 1 0.1 25 nh3 hono 7 newT];
VariableShift(8,:) = [5 3 1 0.1 25 nh3 hono 8 newT];
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
    newT = VariableShift(indexT,9);
  
    serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);
    cd(serial_id_new);
    serial_id3 = sprintf('WaterDynamics_Drying%d_varyingT',desiccationIndex);
    load(strcat(serial_id3,'.mat'));
    waterContentsDynamicsTot{indexT} = waterContentsDynamics;
    waterVolumeDynamicsTot{indexT} = waterVolumeDynamics;
    waterFilmDynamicsTot{indexT} = waterFilmDynamics;
    InvasedIslandDynamicsTot{indexT} = InvasedIslandDynamics;
    gasVolumeDynamicsTot{indexT} = gasVolumeDynamics;  
    serial_idtemp = sprintf('Constant_dark_water_lineardecay_newT%d_NH3%d_HONO%d_Drying%d_varyingT', newT, NH3ppb,HONOppb,desiccationIndex);
    cd(serial_idtemp);
    load('Final.mat','timeLine', 'timeConcDist', 'timeConcDist2', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC','avgWF', 'ActivityCells', 'BioMassDist','patchArea','TotalporeD','porosityM', 'Ymax')
    cd ..
    cd ..
    totSoilVolume = patchArea*TotalporeD./porosityM;
    totSoilWeight = totSoilVolume*1300000;
    AOBActivty = zeros(300,4,1728);
    NOBActivty = zeros(300,4,1728);
    for i = 1:length(timeLine)
        AOBActivty(:,:,i) = ActivityCells{7,i}.*BioMassDist{7,i};
        NOBActivty(:,:,i) = ActivityCells{8,i}.*BioMassDist{8,i};
    end
    %Nitrite production rate by AOB
    NitriteProd = AOBActivty*Ymax(7,7);
    %Nitrite consumption rate by AOB
    NitriteCons = AOBActivty*Ymax(7,8);
    
    totalNitriteReaction = zeros(300,4,1728);
    AmmoniumProd = zeros(300,4,1728);
    NitrateProd = zeros(300,4,1728);
    temp = zeros(300,4);
    temp2 = zeros(300,4);
    AmmoniumConsSum = zeros(1728,1);
    NitrateProdSum = zeros(1728,1);
    NitriteProdSum = zeros(1728,1);
    NitriteConsSum = zeros(1728,1);
    totalpH = zeros(1728,1);
    minpH = zeros(1728,1);
    maxpH = zeros(1728,1);
    meanpH = zeros(1728,1);
    stdpH = zeros(1728,1);
    
    for i = 1:length(timeLine)
        %Total Nitrite reaction
        totalNitriteReaction(:,:,i) = (NitriteProd(:,:,i) + NitriteCons(:,:,i)); % unit g/s per g soil (per patch) local information
        %Ammonium consumption rate by AOB
        AmmoniumProd(:,:,i) = AOBActivty(:,:,i)*Ymax(6,7); % unit g/s per g soil
        %Nitrate production rate by NOB
        NitrateProd(:,:,i) = NOBActivty(:,:,i)*Ymax(5,8); % unit g/s per g soil
        
        temp = AmmoniumProd(:,:,i);
        AmmoniumConsSum(i) = sum(temp(:));
        temp = NitrateProd(:,:,i);
        NitrateProdSum(i) = sum(temp(:));
        temp = NitriteProd(:,:,i);
        NitriteProdSum(i) = sum(temp(:));
        temp = NitriteCons(:,:,i);
        NitriteConsSum(i) = sum(temp(:));
        
        temp = waterVolumeDynamics(:,:,i).*timeConcDist2{4}(:,:,i);
        temp2 = waterVolumeDynamics(:,:,i);
        totalpH(i) = sum(temp(:))./sum(temp2(:));
        minpH(i) = min(min(timeConcDist2{5}(:,:,i)));
        maxpH(i) = max(max(timeConcDist2{5}(:,:,i)));
        meanpH(i) = mean(mean(timeConcDist2{5}(:,:,i)));
        stdpH(i) = std2(timeConcDist2{5}(:,:,i));
        
    end
    
    totVerticalArea = mean(totaldM(:))*Lp*n;
    EffluxAll(:,:,indexT) = effluxList2(:,:)*10^9/dt/plottt./totVerticalArea;
    totalNutriAll(:,:,indexT) = totalNutri(:,:);
    avgWCAll(:,indexT) = avgWC;
    avgWFAll(:,indexT) = avgWF;
    NH4Cons(:,indexT) = AmmoniumConsSum*10^9./totVerticalArea;
    NO3Prod(:,indexT) = NitrateProdSum*10^9./totVerticalArea;
    NO2Cons(:,indexT) = NitriteConsSum*10^9./totVerticalArea;
    NO2Prod(:,indexT) = NitriteProdSum*10^9./totVerticalArea;
    NO2tot(:,indexT) = NO2Cons(:,indexT) + NO2Prod(:,indexT);
    totalpHtot(:,indexT) = totalpH;
    minpHtot(:,indexT) = minpH;
    maxpHtot(:,indexT) = maxpH;
    meanpHtot(:,indexT) = meanpH;
    stdpHtot(:,indexT) = stdpH;
       
end

serial_id2 = sprintf('HT_BSC_photoY3_lineardecay_newT%d_NH3%d_HONO%d_Drying%d', newT, NH3ppb, HONOppb, desiccationIndex);
save(strcat('./Results_HONO_NH3/',serial_id2,'.mat'),'-v7.3');

