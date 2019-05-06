
numberOfSamples = 6;
VariableShift = zeros(numberOfSamples,7);

NH3ppb = 5;
hono = 0.1;
%[0.1 5 10 15]
VariableShift(1,:) = [5 3 1 0.1 25 hono 1];
VariableShift(2,:) = [5 3 1 0.1 25 hono 2];
VariableShift(3,:) = [5 3 1 0.1 25 hono 3];
VariableShift(4,:) = [5 3 1 0.1 25 hono 4];
VariableShift(5,:) = [5 3 1 0.1 25 hono 7];
VariableShift(6,:) = [5 3 1 0.1 25 hono 8];
%VariableShift(8,:) = [5 3 1 0.1 25 3 8];

desiccationIndex = 4;
EffluxAll = zeros(1728,5,numberOfSamples);
totalNutriAll = zeros(1728,8,numberOfSamples);
avgWCAll = zeros(1728,numberOfSamples);
avgWFAll = zeros(1728,numberOfSamples);
alphaCO = 3;

for indexT = 1:numberOfSamples
    examineDay = VariableShift(indexT,1);
    pot1 = VariableShift(indexT,2);
    pCO2frac = VariableShift(indexT,3);
    amountCations = VariableShift(indexT,4);
    avgT = VariableShift(indexT,5);
    HONOppb =VariableShift(indexT,6);
    indexS = VariableShift(indexT,7);

    serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);
    cd(serial_id_new);
    serial_idtemp = sprintf('Constant_dark_water_lineardecay_HONO%d_Drying%d_varyingT', HONOppb, desiccationIndex);
    cd(serial_idtemp);
    load('Final.mat','timeLine', 'totalNutri','totaldM','Lp','n','effluxList2','dt','plottt','avgWC', 'avgWF')
    cd ..
    cd ..
    totVerticalArea = mean(totaldM(:))*Lp*n;
    EffluxAll(:,:,indexT) = effluxList2(:,:)*10^9/dt/plottt./totVerticalArea;
    totalNutriAll(:,:,indexT) = totalNutri(:,:);
    avgWCAll(:,indexT) = avgWC;
    avgWFAll(:,indexT) = avgWF;

end

serial_id2 = sprintf('HT_BSC_photoY3_lineardecay_NH3%d_HONO%d_Drying%d',NH3ppb, HONOppb, desiccationIndex);
save(strcat('./Results_HONO_NH3/',serial_id2,'.mat'),'-v7.3');
