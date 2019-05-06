%function [effluxList2, timeSeriesDeltaM] = EffluxCalculation_dynamics(depthList,timeLine,alphaMatrix,timeConcDist,timeConcDist2,GaseousEffluxMicrobes,InvasedIsland,waterVolume,gasVolume,averageT, amplitudeT, attenuationRate, MaximumIncidence,iniConcentrationGas,sitesC,sitesC2,sitesCg)

[m,n] = size(alphaMatrix);
dataPoints = length(timeLine);

for i = 1:5
    timeSeriesDeltaM{i,1} = zeros(m,n,dataPoints);
end
effluxList2 = zeros(dataPoints,5);
ListGas = [1 2 6 7 8];

for examineT = 1:dataPoints
    if (examineT < preDryHlist + DryHlist+1)&&(examineT > preDryHlist)
        t = (examineT-1)*plottt;
        if rem(t, (60*60/dt)) == 0; %every hour change the matric potential.
            exposedHours = ceil(t*dt/(60*60))+1
        end
        % change conditions
        pot1 = potList(examineT-preDryHlist);
        pot = -pot1/9.81;
        parfor iInd = 1:m
            [WF,percol] = WaterDistAffineHTSatuAProfile(systemAffine(iInd,:,:),pot,Hurst,R1,R2,R);
            waterFilm(iInd,:) = WF(1,:);
            percolProb(iInd,:) = percol(1,:);
        end
        
        [microbeVelocityM, probOccup] = velocityMicrobeMatrix2(waterFilm, pot,V0,R);
        perShare = percolProb*probOccup;
        waterVolume = waterFilm*patchArea; %Volume  of water
        LocalWaterContents = waterFilm./totaldM;
        LocalWaterSaturation = LocalWaterContents./porosityM;
        LocalGasContents = porosityM-LocalWaterContents;
        Aqdiff = (LocalWaterContents.^(10/3)).*porosityM.^(-2);
        localDiffG = 10^4*(LocalGasContents.^(10/3)).*porosityM.^(-2);
        gasVolume = LocalGasContents.*totaldM*patchArea; %Volume of water
        
    end
    
    for i = 1:numberOFN%except sugar and EPS (initial condition for only CO2 and O2)
        if i < 8
            sitesC2{i} = sitesN2{i}./waterVolume;
        end
        sitesC{i} = sitesN{i}./waterVolume;
    end
    sitesC2{5} = -log10(sitesC2{4});
    sitesCa(:,:,1) = sitesCaN1./waterVolume;
    sitesCa(:,:,2) = sitesCaN2./waterVolume;
    
    [intensityListH,temperatureList,MumaxT, HenryConstList{1,1},HenryConstList{2,1},HenryConstList{6,1},HenryConstList{7,1},HenryConstList{8,1}, pKList,Density_air] = EnvironmentProfileDensityAirTmax(depthList,timeLine(examineT),alphaMatrix);
    
    %% Reset Environmental conditions (Light and temperature -> growth rate accordingly)
    thetaTC(:,:) = diag(temperatureList(1,:))*ones(m,n);
    
    pKTotList = zeros(m,n,6); %pK
    for i = 1:6
        pKTotList(:,:,i) = diag(pKList(1,:,i))*ones(m,n);
    end
    
    for i = 3:5
        HenryConstList{i,1} = zeros(1,m);
    end
    rhoAir = zeros(m,n);
    for i = 1:n
        rhoAir(:,i) = Density_air(1,:); %EPS
    end
    
    for iG = 1:length(ListGas)
        gI = ListGas(iG);
        for i = 1:n
            temp2(:,i) = HenryConstList{gI,1}(1,:);
            rhoAir(:,i) = Density_air(1,:); %EPS
        end
        HenryDomain{gI} = temp2;
        AlphaM{gI} = temp2.*waterVolume + gasVolume;
        sitesNGEquil{gI} = PartialPressure(gI)*rhoAir;
        sitesNEquil{gI} = temp2.*sitesCg{gI,1};
    end
    changeN = zeros(m,n,length(timeConcDist));
    
    for i = 1:8
        
        sitesC{i}(:,:) = timeConcDist{i,1}(:,:,examineT);
        sitesCg{i}(:,:) = timeConcDist{i,2}(:,:,examineT);
        changeN(:,:,i) = GaseousEffluxMicrobes{i,examineT}(:,:);
        if i <8
            sitesC2{i}(:,:) = timeConcDist2{i}(:,:,examineT);
        end
        
    end
    
    % InvasedIsland(1,:) = 1;
    %calculate the efflux and equilibrisind the gas-liquid phase: mass
    %conserve
    [efflux, sitesC, sitesCg, sitesC2, deltaNM] = EquilibriumConcentrationFluxDomain(waterVolume,gasVolume,sitesC,sitesC2,sitesCg,HenryDomain,InvasedIsland,AlphaM,sitesNEquil,changeN,iniConcentrationGas);
    effluxList2(examineT,:) = effluxList2(examineT,:) + efflux;% + sum(sum((timeConcDist{1,2}(:,:,examineT)-sitesCg{1,1}.*gasVolume.*InvasedIsland)));
    
    for iD = 1:5
        timeSeriesDeltaM{iD}(:,:,examineT) = deltaNM(:,:,iD);
    end
    
end
%end
