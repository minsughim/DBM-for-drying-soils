numberOfSamples = 9;
VariableShift = zeros(numberOfSamples,9);
examineDay =5;
pot = 3;
examineDay2 = 3;
pCO2frac = 1;
NfixationR = 0.1;
amountCations = 0.1;
avgT = 25;
peneD = 0.2;
reducedLevelC = 0.1;
dryingT = 1;

VariableShift(1,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 1 1];
VariableShift(2,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 1 2];
VariableShift(3,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 1 3];
VariableShift(4,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 0.5 1];
VariableShift(5,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 0.5 2];
VariableShift(6,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 0.5 3];
VariableShift(7,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 0.1 1];
VariableShift(8,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 0.1 2];
VariableShift(9,:) = [examineDay pot pCO2frac NfixationR amountCations avgT peneD 0.1 3];

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
    serial_id = sprintf('HT_BSC_noPG_Tchange_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);
    cd(serial_id);
    %serial_idtemp = sprintf('Constant_light_bubble_Pa_1_2');
    %serial_idtemp = sprintf('Constant_dark_bubble_Pa_1_2');
    serial_idtemp = sprintf('Dinural_cycle_bubble_Pa_1_2');
    % serial_idtemp = sprintf('Dinural_cycle_bubble');

    cd(serial_idtemp);
    load('Final.mat', 'ActivityCells','GaseousEffluxBubbles','LocalGasContents','listV', 'm' ,'n','timeConcDist','timeConcDist2','totaldM','Lp','porosity','dt','plottt','effluxList','popWalkers','Mumean','examineT','timeLine', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesC','sitesC2','sitesCg','popMovie','dormpopWalkers')
    cd ..
    cd ..
    %Efflux Calculation
    [effluxList2, timeSeriesDeltaM] = EffluxCalculationBubble(GaseousEffluxBubbles,depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg);
    
    [effluxList3, timeSeriesDeltaM2] = EffluxCalculation(depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg);
  % serial_id2 = sprintf('HT_BSC_noPG_Tchange_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);

    %save(strcat('./Results/',serial_id,'Light_essence.mat'));
    %save(strcat('./Results/',serial_id,'Dark_essence.mat'));
    %save(strcat('./Results/',serial_id,'Dinural_essence.mat'));
     save(strcat('./Results/',serial_id,'Dinural_Pa_1_2_essence.mat'));
    %save(strcat('./Results/',serial_id,'efflux2.txt'),'effluxList2', '-ascii');
end





