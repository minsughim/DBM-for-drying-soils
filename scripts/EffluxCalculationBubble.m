function [effluxList2, timeSeriesDeltaM] = EffluxCalculationBubble(GaseousEffluxBubbles,depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg)

[m,n] = size(alphaMatrix);
dataPoints = length(timeLine);

for i = 1:5
    timeSeriesDeltaM{i,1} = zeros(m,n,dataPoints);
end
effluxList2 = zeros(dataPoints,5);
ListGas = [1 2 6 7 8];

[intensityList,temperatureList,MumaxT, HenryConstList{1,1},HenryConstList{2,1},HenryConstList{6,1},HenryConstList{7,1},HenryConstList{8,1}, pKList] = EnvironmentProfile_Local(depthList,timeLine,alphaMatrix,averageT, amplitudeT, attenuationRate, MaximumIncidence);
for examineT = 1:dataPoints
    
    for iG = 1:length(ListGas)
        gI = ListGas(iG);
        for i = 1:n
            temp2(:,i) = HenryConstList{gI,1}(examineT,:);
        end
        HenryDomain{gI} = temp2;
        AlphaM{gI} = temp2.*waterVolume + gasVolume;
        sitesNEquil{gI} = temp2.*iniConcentrationGas(gI);
    end
    
    changeN = zeros(m,n,length(timeConcDist));
    
    for i = 1:8
        
        sitesC{i}(:,:) = timeConcDist{i,1}(:,:,examineT);
        sitesCg{i}(:,:) = timeConcDist{i,2}(:,:,examineT);
        changeN(:,:,i) = GaseousEffluxMicrobes{i,examineT}(:,:) + GaseousEffluxBubbles{i,examineT}(:,:);
        if i <8
            sitesC2{i}(:,:) = timeConcDist2{i}(:,:,examineT);
        end
        
    end
    
    % InvasedIsland(1,:) = 1;
    %calculate the efflux and equilibrisind the gas-liquid phase: mass
    %conserve
    [efflux, sitesC, sitesCg, sitesC2, deltaNM] = EquilibriumConcentrationFluxDomain(waterVolume,gasVolume,sitesC,sitesC2,sitesCg,HenryDomain,InvasedIsland,AlphaM,sitesNEquil,changeN,iniConcentrationGas);
    effluxList2(examineT,:) = effluxList2(examineT,:) + efflux;% + sum(sum((timeConcDist{1,2}(:,:,examineT)-sitesCg{1,1}.*gasVolume.*InvasedIsland)));
    
    for iD = 5
        timeSeriesDeltaM{iD}(:,:,examineT) = deltaNM(:,:,iD);
    end
    
end
end
