numberOfSamples = 3;
VariableShift = zeros(numberOfSamples,8);

desiccationIndex = 4;
hono = 0.005;
nh3 = 20;

NdataPoints = 1728;
%[0.1 5 10 15]
VariableShift(1,:) = [5 3 1 0.1 25 nh3 hono 1];
VariableShift(2,:) = [5 3 1 0.1 25 nh3 hono 3];
VariableShift(3,:) = [5 3 1 0.1 25 nh3 hono 5];
%VariableShift(8,:) = [5 3 1 0.1 25 3 8];
EffluxAll = zeros(NdataPoints,5,numberOfSamples);
totalNutriAll = zeros(NdataPoints,8,numberOfSamples);
avgWCAll = zeros(NdataPoints,numberOfSamples);
avgWFAll = zeros(NdataPoints,numberOfSamples);
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
    
    serial_id_new = sprintf('HT_BSC_photoY3_100mu_Day%dPot%.1findex%d',examineDay,pot1,indexS);
    cd(serial_id_new);
     serial_id3 = sprintf('WaterDynamics_Drying%d_varyingT',desiccationIndex);
     load(strcat(serial_id3,'.mat'));
     waterContentsDynamicsTot{indexT} = waterContentsDynamics;
     waterVolumeDynamicsTot{indexT} = waterVolumeDynamics;
     waterFilmDynamicsTot{indexT} = waterFilmDynamics;
     InvasedIslandDynamicsTot{indexT} = InvasedIslandDynamics;
     gasVolumeDynamicsTot{indexT} = gasVolumeDynamics;  
    serial_idtemp = sprintf('Constant_dark_H2ONO_NH3%d_HONO%d_Drying%d_varyingT', NH3ppb,HONOppb,desiccationIndex);
    cd(serial_idtemp);
    load('Final.mat','m','n','timeLine', 'timeConcDist', 'timeConcDist2', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC','avgWF', 'ActivityCells', 'BioMassDist','patchArea','TotalporeD','porosityM', 'Ymax')
    cd ..
    cd ..
    NdataPoints = length(timeLine);
    totSoilVolume = patchArea*TotalporeD./porosityM;
    totSoilWeight = totSoilVolume*1300000;
    AOBActivty = zeros(m,n,NdataPoints);
    NOBActivty = zeros(m,n,NdataPoints);
    for i = 1:length(timeLine)
        AOBActivty(:,:,i) = ActivityCells{7,i}.*BioMassDist{7,i};
        NOBActivty(:,:,i) = ActivityCells{8,i}.*BioMassDist{8,i};
    end
    %Nitrite production rate by AOB
    NitriteProd = AOBActivty*Ymax(7,7);
    %Nitrite consumption rate by AOB
    NitriteCons = AOBActivty*Ymax(7,8);
    
    totalNitriteReaction = zeros(m,n,NdataPoints);
    AmmoniumProd = zeros(m,n,NdataPoints);
    NitrateProd = zeros(m,n,NdataPoints);
    temp = zeros(m,n);
    temp2 = zeros(m,n);
    AmmoniumConsSum = zeros(NdataPoints,1);
    NitrateProdSum = zeros(NdataPoints,1);
    NitriteProdSum = zeros(NdataPoints,1);
    NitriteConsSum = zeros(NdataPoints,1);
    totalpH = zeros(NdataPoints,1);
    minpH = zeros(NdataPoints,1);
    maxpH = zeros(NdataPoints,1);
    meanpH = zeros(NdataPoints,1);
    
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
       
end

serial_id2 = sprintf('HT_BSC_photoY3_100mu_H2ONO_NH3%d_HONO%d_Drying%d',NH3ppb, HONOppb, desiccationIndex);
save(strcat('./Results_HONO_NH3/',serial_id2,'.mat'),'-v7.3');

