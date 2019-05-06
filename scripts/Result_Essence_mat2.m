numberOfSamples = 3;
VariableShift = zeros(numberOfSamples,9);
examineDay =3;
pot1 = 3;
examineDay2 = 5;
pCO2frac = 1;
NfixationR = 0.1;
amountCations = 0.1;
avgT = 25;
peneD = 0.2;
alphaCO = 5;
reducedLevelC = 0.1;
dryingT = 1;

for indexT = 1:numberOfSamples
    
    indexS = indexT;
    serial_id = sprintf('HT_BSC_photoY2_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);
    cd(serial_id);
    %serial_idtemp = sprintf('Constant_light');
    %serial_idtemp = sprintf('Constant_dark');
    serial_idtemp = sprintf('Dinural_cycle');
    % serial_idtemp = sprintf('Dinural_cycle_bubble');
    cd(serial_idtemp);
    %load('Final.mat', 'BioMassDist', 'ActivityCells','GaseousEffluxBubbles','LocalGasContents','listV', 'm' ,'n','timeConcDist','timeConcDist2','totaldM','Lp','porosity','dt','plottt','effluxList','popWalkers','Mumean','examineT','timeLine', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesC','sitesC2','sitesCg','popMovie','dormpopWalkers')
    load('Final.mat', 'BioMassDist', 'ActivityCells','LocalGasContents','listV', 'm' ,'n','timeConcDist','timeConcDist2','totaldM','Lp','porosity','dt','plottt','effluxList','popWalkers','Mumean','examineT','timeLine', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesC','sitesC2','sitesCg','popMovie','dormpopWalkers')
    cd ..
    cd ..
    %Efflux Calculation
%    [effluxList2, timeSeriesDeltaM] = EffluxCalculationBubble(GaseousEffluxBubbles,depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg);    
    [effluxList2, timeSeriesDeltaM] = EffluxCalculation(depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg);
  % serial_id2 = sprintf('HT_BSC_noPG_Tchange_Day%dPot%.1fpCO2frac%dNfix%.1fCat%.1favT%dPeneD%.1fredC%.2findex%d',examineDay,pot1,pCO2frac,NfixationR, amountCations, avgT, peneD, reducedLevelC, indexS);

    %save(strcat('./Results/',serial_id,'Light_essence.mat'));
    %save(strcat('./Results/',serial_id,'Dark_essence.mat'));
    save(strcat('./Results2/',serial_id,'_Dinural_essence.mat'));
  %   save(strcat('./Results2/',serial_id,'essence.mat'));
    %save(strcat('./Results/',serial_id,'efflux2.txt'),'effluxList2', '-ascii');
end





