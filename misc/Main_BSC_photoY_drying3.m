function Main_BSC_photoY_drying3(lightOn2,examineDay,dryingT, pCO2frac, amountCations, avgT, alphaCO, pot1, plottt, indexS,maximumThreadt)
rng('shuffle')

%examineDay =5;
%dryingT = 5;
%NfixationR = 0.1;
%pCO2frac = 1;
%amountCations = 0.1;
%avgT = 25;
%peneD = 0.2;
%pot1 = 3;
%plottt = 5;
%indexS = 9;
%alphaCO = 1;
%maximumThreadt = 2;
%lightOn = 0;

global m n Rad3 Rad3Over2 dt Lp dl patchArea  R1 R2 %modelling parameters
global V0 chi0 mumax Ks Ymax mrate Vu R rho % prameters for microbes
global dataPoints ed listE numberOFN numberOFP testMinT Hurst
global averageT amplitudeT attenuationRate MaximumIncidence repsiRatio
global iniConcentrationGas iniConcentration

serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);
cd(serial_id_new);
serial_idtemp_dd = sprintf('Dinural_cycle');
cd(serial_idtemp_dd);
load('Final.mat')
cd ..
if lightOn2 == 1
    serial_idtemp = sprintf('Constant_light_Drying%ddays_varyingT', dryingT);
    amplitudeT = 0; %experiemnt with the constnat condition (with average velocity and maximum incidence of light)
else if lightOn2 == 0
        serial_idtemp = sprintf('Constant_dark_Drying%ddays_varyingT', dryingT);
        amplitudeT = 0; %experiemnt with the constnat condition (with average velocity and maximum incidence of light)
    else
        serial_idtemp = sprintf('Dinural_cycle_Drying%ddays_varyingT', dryingT);
    end
end
mkdir(serial_idtemp);
cd ..
%% Changes in water requiredment for phototrophs
Parameters_Modification_water_shortage

%% Drying
examineDay = dryingT;
reducedLevelC = 10^(-16); %to start with zero input from bottom, it increases slowly depending on the hydration condition
satuIndex = 0; % saturation index 1: indicating completely saturated case, saturation index 0: indicating field capacit corresponds to 0.66 saturation condition
repsiRatio = 0.1; %
Zratio = amountCations; %nut pH = 7.5;
Zion = Zratio*100*10^(-3)/40; %in mol/L (net ion concentration that t
calciumF = amountCations;
initCalcium = calciumF*100*10^(-3)/40; %%100 maicro g/g

%Change EPS production
%EPSfraction = 0.1;

%% Information for bubble formation
pN2endOthers = 1 - sum(PartialPressure(1:8)); %non-reactive gaseous elements
totNgBubble = zeros(m,n,numberOFN);
P_a = 1.4;

%% Time information for the fixed time step
dt = 60; %second %1 min.
testMinT = 0.1; % second: diffusion equation resolution of 1 second.

examineH = dryingT*24 + 24*2; %sustain as the dry soil for two days
timeLine = (dt*plottt)/3600:(dt*plottt)/3600:examineH;
T = (3600/dt)*examineH; %examine time, (seconds)*hours
Td = (3600/dt)*dryingT*24; %examine time, (seconds)*hours
dataPoints = T/plottt;
timeLapseList = zeros(dataPoints,1);
DeltaT = dt;
DryHlist = Td/plottt;
potList = linspace(0.5,30,dryingT*24);
%% Preallocate variables
endIndex = length(Mumean);
clear Mumean

popWalkers = zeros(dataPoints,numberOFP);
dormpopWalkers = zeros(dataPoints,numberOFP);
totalBioMass = zeros(dataPoints,numberOFP);

ActivityCells = cell(numberOFP,dataPoints);
Mumean= zeros(dataPoints,numberOFP);
popMovie= cell(numberOFP,dataPoints);
dormpopMovie= cell(numberOFP,dataPoints);
measureMMovie= cell(numberOFP,dataPoints);
totalNutriPart = cell(numberOFN,1);
tempsitesC = cell(numberOFN,1);
changeNurient = cell(numberOFN,1);
sitesN = cell(numberOFN,1); % mass of nutrients
sitesN2 = cell(numberOFN,1); % mass of nutrients
%sitesCg = cell(numberOFN,1); %concentration of nutrients
sitesNg = cell(numberOFN,1); %concentration of nutrients

reactionC = zeros(m,n,numberOFN+3);
totalNutri = zeros(dataPoints,numberOFN);

for i = 1:8
    sitesN{i} = sitesC{i}.*waterVolume;
end

for i = 1:7
    sitesN2{i} = sitesC2{i}.*waterVolume;
end
%ZionMN = ZionM.*waterVolume;
sitesCaN1 = sitesCa(:,:,1).*waterVolume;
sitesCaN2 = sitesCa(:,:,2).*waterVolume;

timeConcDist = cell(numberOFN,1);
timeConcDist2 = cell(9,1); %EPS, CO3, NH3, H+, PH

for i = 1:numberOFN%except sugar and EPS (initial condition for only CO2 and O2)
    changeNurient{i,1} = zeros(m+2, n+2);
    totalNutriPart{i,1} = zeros(plottt,1);
    timeConcDist{i,1} = zeros(m, n,dataPoints);
    timeConcDist{i,2} = zeros(m, n,dataPoints);
    if i < 8
        timeConcDist2{i,1} = zeros(m, n,dataPoints);
        sitesC2{i} = sitesN2{i}./waterVolume;
    end
    sitesC{i} = sitesN{i}./waterVolume;
end

sitesC2{5} = -log10(sitesC2{4});
%ZionM = ZionMN./waterVolume;

effluxList = zeros(dataPoints,5);
ListGas = [1 2 6 7 8];
AlphaM = cell(numberOFN,1);
sitesNEquil = cell(numberOFN,1);
sitesNGEquil = cell(numberOFN,1);
HenryDomain = cell(numberOFN,1);
temp2 = zeros(m,n);
for iG = 1:length(ListGas)
    gI = ListGas(iG);
    for i = 1:n
        temp2(:,i) = HenryConstList{gI,1}(1,:);
    end
    HenryDomain{gI} = temp2;
    AlphaM{gI} = temp2.*waterVolume + gasVolume;
    sitesNGEquil{gI} = sitesCg{gI,1};
    sitesNEquil{gI} = temp2.*sitesCg{gI,1};
end

tic

day = 1
%% Starting dynamics
timeM = zeros(m,n);

for examineT = 1:dataPoints
    if examineT < DryHlist+1
        
        t = (examineT-1)*plottt;
        if rem(t, (60*60/dt)) == 0; %every hour change the matric potential.
            exposedHours = ceil(t*dt/(60*60))+1
            % change conditions
            pot1 = potList(exposedHours);
            pot = -pot1/9.81;
            parfor iInd = 1:m
                systemAffineFrac = zeros(1,n,2);
                [WF,percol,WS,CP,SintA,LengthL] = WaterDistAffineHTSatuAProfile(systemAffine(iInd,:,:),pot,Hurst,R1,R2,R);
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
            
            % new boundary conditions
            [realIsland,town2,results2,numberOfCluster] = IslandStatHex(LocalGasContents,0.2); %Invasive percolation for BCC case Pc = 0.2
            
            EpsilonG = mean(LocalGasContents(:));
            Phi = mean(porosityM(:));
            tau = log((2*EpsilonG*EpsilonG+0.04*EpsilonG)/Phi/Phi)/log(EpsilonG/Phi);
            ThetaV = mean(LocalWaterContents(:));
            denom = (Phi*Phi*(EpsilonG/Phi)^tau)/(EpsilonG+6.42e-3*ThetaV);
            BreakThroughGas = zeros(m,n);
            for i = 1:length(depthList)
                loc = depthList(i)*1000;
                timeM = timeM+town2*plottt*dt;
                for j = 1:n                   
                    if timeM(i,j)==0;
                        BreakThroughGas(i,j) = 0;
                    else       
                        BreakThroughGas(i,j) = erfc(loc./(0.23*denom*timeM(i,j)));
                    end
                end
            end
            
            InvasedIsland = town2.*BreakThroughGas;
            for i = 1:numberOfCluster
                ListLocation = find(realIsland==i);
                InvasedIsland(ListLocation) = (gasVolume(ListLocation)'*InvasedIsland(ListLocation))/sum(gasVolume(ListLocation));
            end
            InvasedIsland(1,:) = 1;
            
            ListUnsatu = find(waterSatu<1);
            if isempty(ListUnsatu)==1
                exRatewithoutDiff =zeros(m,n);
            else
                exRatewithoutDiff = ((localDiffG./(1-LocalWaterSaturation))-(Aqdiff./LocalWaterSaturation)).*sepecifInterA./TotalporeD;
            end
            
            lambdaMatrix = ((1-porosityM)/2.9 + LocalWaterContents/0.57 + LocalGasContents/0.0012).^(-1);%Soil Thermal conductivity [W/mK]
            CvMatrix = 10^(6)*(1.94*(1-porosityM)+4.189*LocalWaterContents+ 0.0012*LocalGasContents); %Soil Thermal diffusivity[m^2/s]
            alphaMatrix = lambdaMatrix./CvMatrix; %Soil Thermal diffusivity[m^2/s]
            
        end
        
    end
    
    
    for i = 1:numberOFN%except sugar and EPS (initial condition for only CO2 and O2)
        if i < 8
            sitesC2{i} = sitesN2{i}./waterVolume;
        end
        sitesC{i} = sitesN{i}./waterVolume;
    end
    sitesC2{5} = -log10(sitesC2{4});
    %ZionM = ZionMN./waterVolume;
    sitesCa(:,:,1) = sitesCaN1./waterVolume;
    sitesCa(:,:,2) = sitesCaN2./waterVolume;
    
    [intensityListH,temperatureList,MumaxT, HenryConstList{1,1},HenryConstList{2,1},HenryConstList{6,1},HenryConstList{7,1},HenryConstList{8,1}, pKList,Density_air] = EnvironmentProfileDensityAirTmax(depthList,timeLine(examineT),alphaMatrix);
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
    %% water informaion update
    tempwaterVolume = generateBCcylinderY(waterVolume);
    tempWFm = tempwaterVolume/patchArea;
    tempWFm(1:2,:) = tempWFm(1:2,:)./100; %top of the surface has lower effecive water film thickness
    normWater = cell(m,n);
    for y = 1:m
        for x = 1:n
            normWater{y,x} = WeightFactorWperiodCylinder(tempWFm(y:y+2,x:x+2),x,y);
        end
    end
    for i = 1:numberOFN
        totalNutriPart{i,1} = zeros(plottt,1);
    end
    
    %% Reset Environmental conditions (Light and temperature -> growth rate accordingly)
    intensityProfile = diag(intensityList(examineT,:))*ones(m,n);
    MumaxTt(:,:) = diag(MumaxT(1,:))*ones(m,n);
    thetaTC(:,:) = diag(temperatureList(1,:))*ones(m,n);
    
    pKTotList = zeros(m,n,6); %pK
    for i = 1:6
        pKTotList(:,:,i) = diag(pKList(1,:,i))*ones(m,n);
    end
    
    %% Population update
    PopulationMapS = cell(m,n);
    PopulationMovie = cell(numberOFP,plottt);
    dormPopulationMapS = cell(m,n);
    dormPopulationMovie = cell(numberOFP,plottt);
    muList = cell(numberOFP,1);
    PopulationPart = cell(numberOFP,1);
    dormPopulationPart = cell(numberOFP,1);
    for i = 1:numberOFP
        for j=1:plottt
            PopulationMovie{i,j} = sparse(m,n);
            dormPopulationMovie{i,j} = sparse(m,n);
        end
        muList{i} = zeros(m,n,plottt);
    end
    popWalkersPart = zeros(plottt,numberOFP);
    dormpopWalkersPart = zeros(plottt,numberOFP);
    MumeanPart = zeros(plottt,numberOFP);
    wakeUpTemp = ones(m,n);
    
    totCsum = 0*totCsum;
    totNgsum = 0*totNgsum;
    changeN = 0*changeN ;
    reactionC = 0*reactionC;
    totNgBubble = 0*totNgBubble;
    
    
    for tt = 1:plottt
        
        t = (examineT-1)*plottt + tt;
        
        if rem(t, (6*60*60/dt)) == 0; %every 6 hours save
            exposedHours = t*dt/(60*60)
            serial_id3 = sprintf('HTBioCrustEPSCNcycleHour%d',exposedHours);
            save(strcat('./',serial_id_new,'/',serial_idtemp,'/',serial_id3,'.mat'),'porosity','totaldM','LocalGasContents','listV','BioMassDist','GaseousEffluxBubbles','ActivityCells','popWalkers', 'dormpopWalkers', 'popMovie','dormpopMovie','totalBioMass','timeConcDist','timeConcDist2','effluxList','timeLine','sitesC','sitesC2', 'sitesCa','totNutConsumption', 'GaseousEffluxMicrobes', 'Mumean','dataPoints', 'totalNutri','exposedHours','Lp','m','n','numberOFP','numberOFN','averageT','depthList','dt','plottt','examineT', 'depthList','alphaMatrix','GaseousEffluxMicrobes','InvasedIsland','waterVolume','gasVolume','averageT', 'amplitudeT', 'attenuationRate', 'MaximumIncidence','iniConcentrationGas','sitesCg');
        end
        
        for i = 1:numberOFN
            sitesN{i,1} = waterVolume.*sitesC{i,1};
            if i < 8
                sitesN2{i,1} = waterVolume.*sitesC2{i,1};
            end
        end
        
        %calculate the efflux and equilibrisind the gas-liquid phase: mass
        %conserve
        [efflux, sitesC, sitesCg, sitesC2] = EquilibriumConcentrationDensity(waterVolume,gasVolume,sitesC,sitesC2,sitesCg,HenryDomain,InvasedIsland,AlphaM,sitesNEquil,sitesNGEquil,changeN);
        effluxList(examineT,:) = effluxList(examineT,:) + efflux;% + sum(sum((timeConcDist{1,2}(:,:,examineT)-sitesCg{1,1}.*gasVolume.*InvasedIsland)));
        
        gradApparent = cell(numberOFP,2);
        possibleDist = cell(numberOFP,2);
        jumpTProbM = cell(numberOFP,2);
        wakeUp = 0*wakeUp;
        
        %% calcualting mumax based on the new environmental conditions
        [apparentGrowth(:,:,1:4),respiGrowthR] = PhotoGrowth(sitesC, sitesC2, MumaxTt,intensityProfile);
        [apparentGrowth(:,:,5),apparentGrowth(:,:,6),apparentGrowth(:,:,7),apparentGrowth(:,:,8)] = ChemiGrowthCrust(sitesC,sitesC2,MumaxTt);
        
        for i = 1:numberOFP
            wakeUp(:,:,i) = wakeUpTemp.*(apparentGrowth(:,:,i)>mrate(i));
        end
        
        %% update walkers infromation
        clear PopulationMapS dormPopulationMapS
        PopulationMapS = cell(m,n);
        dormPopulationMapS = cell(m,n);
        maxNindiv = size(walkerHist,2);
        walkerHistTemp = cell(1,maxNindiv);
        nIndv = 0;
        
        PopulationS = 0*PopulationS;
        
        dormPopulationS = 0*dormPopulationS;
        testV = 0*testV;
        testVinactive = 0*testVinactive;
        measVtot = 0*measVtot;
        pop = 0*pop;
        
        for i = 1:maxNindiv
            testWalker = walkerHist{i};
            switch testWalker.status
                case 1
                    nIndv = nIndv +1;
                    position =testWalker.positionS;
                    xx = position(1);
                    yy = position(2);
                    speciesNum = testWalker.sp;
                    PopulationMovie{speciesNum,tt}(xx,yy) = PopulationMovie{speciesNum,tt}(xx,yy)+1;
                    PopulationS(xx,yy) = PopulationS(xx,yy)+ 1;
                    PopulationMapS{xx,yy}(PopulationS(xx,yy)) = nIndv;
                    testWalker.number = nIndv;
                    walkerHistTemp{nIndv} = testWalker;
                    testV(xx,yy,speciesNum) = testV(xx,yy,speciesNum) + walkerHist{i}.V; %Total volume of cells at (xx,yy) for each species
                    pop(xx,yy,speciesNum) = pop(xx,yy,speciesNum) + 1;  %Total populations of cells at (xx,yy) for each species
                case 2
                    nIndv = nIndv +1;
                    position =testWalker.positionS;
                    xx = position(1);
                    yy = position(2);
                    speciesNum = testWalker.sp;
                    Wake = wakeUp(xx,yy,speciesNum)*(exp(-(testWalker.waitingT-(germT*6*60*dt))/(germT*60*60*dt))<rand);
                    if Wake == 1
                        %if  wakeUp{speciesNum,1}(xx,yy) == 1 %No germination time just become active immidiately
                        PopulationMovie{speciesNum,tt}(xx,yy) = PopulationMovie{speciesNum,tt}(xx,yy)+1;
                        PopulationS(xx,yy) = PopulationS(xx,yy)+ 1;
                        PopulationMapS{xx,yy}(PopulationS(xx,yy)) = nIndv;
                        testWalker.number = nIndv;
                        testWalker.status = 1;
                        testWalker.velocity = microbeVelocityM(xx,yy);
                        walkerHistTemp{nIndv} = testWalker;
                        testV(xx,yy,speciesNum) = testV(xx,yy,speciesNum) + walkerHist{i}.V; %Total volume of cells at (xx,yy) for each species
                        pop(xx,yy,speciesNum) = pop(xx,yy,speciesNum) + 1;
                    else
                        dormPopulationMovie{speciesNum,tt}(xx,yy) = dormPopulationMovie{speciesNum,tt}(xx,yy)+1;
                        dormPopulationS(xx,yy) = dormPopulationS(xx,yy)+ 1;
                        dormPopulationMapS{xx,yy}(dormPopulationS(xx,yy)) = nIndv;
                        testWalker.number = nIndv;
                        testWalker.waitingT = testWalker.waitingT + dt*(wakeUp(xx,yy,speciesNum)==1);
                        walkerHistTemp{nIndv} = testWalker;
                        testVinactive(xx,yy,speciesNum) = testVinactive(xx,yy,speciesNum) + walkerHist{i}.V; %Total volume of cells at (xx,yy) for each specie
                        measVtot(xx,yy,speciesNum) = measVtot(xx,yy,speciesNum) + rho*walkerHist{i}.V*mrate(speciesNum)*exp(-testWalker.waitingT*mrate(speciesNum));
                    end
            end
        end
        
        OccuationDensity = (PopulationS+dormPopulationS)./MaxCapacityPop;
        
        HighDensityPatch = (OccuationDensity>1).*OccuationDensity + (OccuationDensity<=1).*ones(m,n);
        numberOfWalkers = nIndv;
        clear walkerHist;
        walkerHist = walkerHistTemp(1:numberOfWalkers);
        clear walkerHistTemp;
        
        %% Calculate local consumption and expected growth rate for each species
        reactionC = 0*reactionC;
        
        tempCN(:,:,1) = ones(m,n) - YmaxCNratio(1).*sitesC{5}./sitesC{2};
        tempCN(:,:,2) = ones(m,n) - YmaxCNratio(2).*sitesC{6}./sitesC{2};
        tempCN(:,:,3) = ones(m,n) - YmaxCNratio(3).*sitesC{5}./sitesC{3};
        tempCN(:,:,4) = ones(m,n) - YmaxCNratio(4).*sitesC{6}./sitesC{3};
        CNratioX = NfixationR*(1+0.5*tanh(tempCN));
        
        totalGrowthMass = apparentGrowth.*testV*rho; %photosynthesis + others
        for y = 1:m
            for x = 1:n
                if isempty(PopulationMapS{y,x}) == 0
                    for i = 1:4
                        YmaxCandN(:,i) = (1-CNratioX(y,x,i))*Ymax(:,i)+CNratioX(y,x,i)*YN2fix';
                    end
                    temp = zeros(numberOFP,1);
                    temp(:) = totalGrowthMass(y,x,:);
                    reactionC(y,x,:) = YmaxCandN*temp;
                else
                    reactionC(y,x,:) = 0;
                end
            end
        end
        
        %% Water might be not enough for the photosynthesis: Then shutdown photosynthesis
        
        reactionCwater = abs(reactionC(:,:,10)).*(reactionC(:,:,10)<0);
        WaterShortage = (reactionCwater<waterVolume*10^(3)).*(1-reactionCwater./waterVolume*10^(-3)); %change in water volume larger than 0.1% than, consider it as water shortage.
        
        % reevaluate photogrowth term based on water shortage
        for i = 1:4
            apparentGrowth(:,:,i) = WaterShortage.*apparentGrowth(:,:,i);
        end
        % Nitrite oxidiser also requires water
        apparentGrowth(:,:,8) = WaterShortage.*apparentGrowth(:,:,8);
        
        
        totalGrowthMass = apparentGrowth.*testV*rho;
        RespiGrowthMass = respiGrowthR.*testV*rho;
        
        
        % Water contents dependent EPS production rate:: EPS fraction
        HydroEPSM = EPSfraction./(1+exp((sitesC2{1}-EPSCcrit)./(EPSCcrit*LocalWaterContents)));
        
        for y = 1:m
            for x = 1:n
                if isempty(PopulationMapS{y,x}) == 0
                    for i = 1:4
                        YmaxCandN(:,i) = (1-CNratioX(y,x,i))*Ymax(:,i)+CNratioX(y,x,i)*YN2fix';
                    end
                    temp = zeros(numberOFP,1);
                    temp(:) = totalGrowthMass(y,x,:);
                    reactionC(y,x,:) = YmaxCandN*temp;
                    reactionCEPS(y,x) = HydroEPSM(y,x)*YmaxCandN(4,1:4)*temp(1:4);
                    temp(:) = RespiGrowthMass(y,x,:);
                    reactionC(y,x,:) = reactionC(y,x,:) + reshape(YrespiTot*temp,[1,1,11]);
                else
                    reactionC(y,x,:) = 0;
                end
            end
        end
        
        reactionC(:,:,4) = reactionC(:,:,4)-reactionCEPS;
        
        diffMatrix = 0*diffMatrix;
        diffMatrix(:,:,1) = DiffusionList(1).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,2) = DiffusionList(2).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,3) = DiffusionList(3).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,4) = DiffusionList(4).*Aqdiff.*max(exp(-0.05*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,5) = DiffusionList(5).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,6) = DiffusionList(6).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,7) = DiffusionList(7).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix(:,:,8) = DiffusionList(8).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        
        parfor i = 1:8
            switch i
                case 1 %O2
                    Diffresult{i} = DiffusionProcess_aq_gas_Invase(diffMatrix(:,:,i),waterFilm,waterVolume,gasVolume,sitesC{i}(:,:),sitesCg{i}(:,:),HenryDomain{i}(:,:),reactionC(:,:,i),InvasedIsland,sitesNEquil{i}(:,:),sitesNGEquil{i}(:,:),AlphaM{i}(:,:),listE,Lp,dt,testMinT);
                case 2 %CO2
                    Diffresult{i} = DiffusionProcess_aq_gas_Invase(diffMatrix(:,:,i),waterFilm,waterVolume,gasVolume,sitesC{i}(:,:),sitesCg{i}(:,:),HenryDomain{i}(:,:),reactionC(:,:,i),InvasedIsland,sitesNEquil{i}(:,:),sitesNGEquil{i}(:,:),AlphaM{i}(:,:),listE,Lp,dt,testMinT);
                case 3 %HCO3-
                    Diffresult{i} = DiffusionProcess_aq(diffMatrix(:,:,i),waterFilm,waterVolume,sitesC{i}(:,:),reactionC(:,:,i),listE,Lp,dt,testMinT);
                case 4 %CH2O
                    Diffresult{i} = DiffusionProcess_aq(diffMatrix(:,:,i),waterFilm,waterVolume,sitesC{i}(:,:),reactionC(:,:,i),listE,Lp,dt,testMinT);
                case 5 %NO3-
                    Diffresult{i} = DiffusionProcess_aq(diffMatrix(:,:,i),waterFilm,waterVolume,sitesC{i}(:,:),reactionC(:,:,i),listE,Lp,dt,testMinT);
                case 6 %NH4+
                    Diffresult{i} = DiffusionProcess_aq(diffMatrix(:,:,i),waterFilm,waterVolume,sitesC{i}(:,:),reactionC(:,:,i),listE,Lp,dt,testMinT);
                case 7 %NO2-
                    Diffresult{i} = DiffusionProcess_aq(diffMatrix(:,:,i),waterFilm,waterVolume,sitesC{i}(:,:),reactionC(:,:,i),listE,Lp,dt,testMinT);
                case 8 %N2O
                    Diffresult{i} = DiffusionProcess_aq_gas_Invase(diffMatrix(:,:,i),waterFilm,waterVolume,gasVolume,sitesC{i}(:,:),sitesCg{i}(:,:),HenryDomain{i}(:,:),reactionC(:,:,i),InvasedIsland,sitesNEquil{i}(:,:),sitesNGEquil{i}(:,:),AlphaM{i}(:,:),listE,Lp,dt,testMinT);
            end
        end
        
        diffMatrix2(:,:,1) = DiffusionAuxil(1).*Aqdiff.*max(exp(-1*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix2(:,:,2) = DiffusionAuxil(2).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix2(:,:,3) = DiffusionAuxil(3).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix2(:,:,4) = DiffusionAuxil(6).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        diffMatrix2(:,:,5) = DiffusionAuxil(7).*Aqdiff.*max(exp(-0.02*(sitesC2{1}/1000).^(0.5)),10^(-4));
        
        
        parfor i = 1:5
            switch i
                case 1 %EPS
                    Diffresult2{i} = DiffusionProcess_aq(diffMatrix2(:,:,i),waterFilm,waterVolume,sitesC2{i}(:,:),reactionCEPS,listE,Lp,dt,testMinT);
                case 2 %CO32-
                    Diffresult2{i} = DiffusionProcess_aq(diffMatrix2(:,:,i),waterFilm,waterVolume,sitesC2{i}(:,:),zeros(m,n),listE,Lp,dt,testMinT);
                case 3 %NH3
                    Diffresult2{i} = DiffusionProcess_aq_gas_Invase(diffMatrix2(:,:,i),waterFilm,waterVolume,gasVolume,sitesC2{i}(:,:),sitesCg{6}(:,:),HenryDomain{6}(:,:),zeros(m,n),InvasedIsland,sitesNEquil{6}(:,:),sitesNGEquil{6}(:,:),AlphaM{6}(:,:),listE,Lp,dt,testMinT);
                case 4 %HONO
                    Diffresult2{i} = DiffusionProcess_aq_gas_Invase(diffMatrix2(:,:,i),waterFilm,waterVolume,gasVolume,sitesC2{6}(:,:),sitesCg{7}(:,:),HenryDomain{7}(:,:),zeros(m,n),InvasedIsland,sitesNEquil{7}(:,:),sitesNGEquil{7}(:,:),AlphaM{7}(:,:),listE,Lp,dt,testMinT);
                case 5 %Ca
                    Diffresult2{i} = DiffusionProcess_aq(diffMatrix2(:,:,i),waterFilm,waterVolume,sitesC2{7}(:,:),zeros(m,n),listE,Lp,dt,testMinT);
            end
        end
        
        for i = 1:numberOFN
            totCsum(:,:,i) = totCsum(:,:,i) + Diffresult{i}.totConsumed;
            totNgsum(:,:,i) = totNgsum(:,:,i) + Diffresult{i}.changeN;
            sitesC{i} = Diffresult{i}.sitesC;
            sitesCg{i} = Diffresult{i}.sitesCg;
            nutShare(:,:,i) = Diffresult{i}.nutShare;
        end
        
        sitesC2{1} = Diffresult2{1}.sitesC;
        sitesC2{2} = Diffresult2{2}.sitesC;
        sitesC2{3} = Diffresult2{3}.sitesC;
        sitesCg{6} = Diffresult2{3}.sitesCg;
        sitesC2{6} = Diffresult2{4}.sitesC;
        sitesCg{7} = Diffresult2{4}.sitesCg;
        sitesC2{7} = Diffresult2{5}.sitesC;
        
        
        %% Bubbles escape when the built up partial pressure of gasous elements are larger than the fraction P_a
        
        temp = pN2endOthers*ones(m,n);
        for iG = 1:length(ListGas)
            gI = ListGas(iG);
            if gI == 6
                temp = temp + sitesC2{3}./HenryDomain{gI}./rhoAir; %Partial pressure divided by the atmpspheric pressure
            else if gI == 7
                    temp = temp + sitesC2{4}./HenryDomain{gI}./rhoAir;
                else
                    temp = temp + sitesC{gI}./HenryDomain{gI}./rhoAir;  %temp = total partial pressure including all the gaseous elements
                end
            end
        end
        
        for iG = 1:length(ListGas)
            gI = ListGas(iG);
            if gI == 6
                sitesCtemp = (temp > P_a).*sitesNEquil{gI} + (temp <= P_a).*sitesC2{3};
                totNgBubble(:,:,gI) = totNgBubble(:,:,gI) + (sitesC2{3} - sitesCtemp).*waterVolume;
                sitesC2{3} = sitesCtemp;
            else if gI == 7
                    sitesCtemp = (temp > P_a).*sitesNEquil{gI} + (temp <= P_a).*sitesC2{4};
                    totNgBubble(:,:,gI) = totNgBubble(:,:,gI) + (sitesC2{4} - sitesCtemp).*waterVolume;
                    sitesC2{4} = sitesCtemp;
                else
                    sitesCtemp = (temp > P_a).*sitesNEquil{gI} + (temp <= P_a).*sitesC{gI};
                    totNgBubble(:,:,gI) = totNgBubble(:,:,gI) + (sitesC{gI} - sitesCtemp).*waterVolume;
                    sitesC{gI} = sitesCtemp;
                end
            end
            
        end
        
        %% calcualting new growth rate based on nutshare: nutrient limitation
        if min(nutShare(:))~=1
            %limitF1 = min(nutShare(:))
            %limitF2 = max(nutShare(:))
            %DeltaT
            %Effect of pH applies to the next time step: Here only mass
            %share is applied
            [apparentGrowth(:,:,1:4),respiGrowthR] = PhotoGrowthShare(sitesC, sitesC2,MumaxTt,intensityProfile,nutShare);
            [apparentGrowth(:,:,5),apparentGrowth(:,:,6),apparentGrowth(:,:,7),apparentGrowth(:,:,8)] = ChemiGrowthCrustShare(sitesC, sitesC2, MumaxTt,nutShare);
        end
        
        % reevaluate photogrowth term based on water shortage
        for i = 1:4
            apparentGrowth(:,:,i) = WaterShortage.*apparentGrowth(:,:,i);
        end
        
        % Nitrite oxidiser also requires water
        apparentGrowth(:,:,8) = WaterShortage.*apparentGrowth(:,:,8);
        totalGrowth = apparentGrowth + respiGrowthR;
        
        %% Update possible displacement and ju mping probability based on the grwoth rate field
        for i = 1:numberOFP
            gradApparent{i,1} = CalculateGradientField(totalGrowth(:,:,i),mumax(i),m,n,Lp);
            [possibleDist{i,1}, jumpTProbM{i,1}] = MicrobeExpMobParBioCrustSpaceConstrain(HighDensityPatch,(PopulationMovie0(:,:,i)>0).*microbeVelocityM, normWater, gradApparent{i,1}, perShare);
        end
        
        %% Update individuals based on spatial information
        
        muSp = 0*muSp;
        %% Update individuals based on spatial information
        muSp = zeros(m,n,numberOFP);
        for y = 1:m
            for x = 1:n
                if isempty(PopulationMapS{y,x}) == 0
                    for i =1:numberOFP
                        muSp(y,x,i) = totalGrowth(y,x,i) - mrate(i);
                        muList{i}(y,x,tt) = totalGrowth(y,x,i)*pop(y,x,i);
                    end
                    listOfWalkers = PopulationMapS{y,x}(:);
                    for ii = 1:length(listOfWalkers)
                        index = listOfWalkers(ii);
                        speciesNum = walkerHist{index}.sp;
                        walkerHist{index}.muGcorr = muSp(y,x,speciesNum);
                        walkerHist{index}.positionV = walkerHist{index}.positionV + possibleDist{speciesNum,1}{y,x}*DeltaT;
                        walkerHist{index}.fluxProb = walkerHist{index}.fluxProb + jumpTProbM{speciesNum,1}{y,x};
                    end
                end
            end
        end
        
        %% Matrix for Parallel caclaution
        for i = (numberOfWalkers+1):(maximumThreadt*ceil(numberOfWalkers/maximumThreadt))
            walkerHist{i}.status = 0;
        end
        
        walkerHistM = reshape(walkerHist,[],maximumThreadt);
        ListNewWalkersM = zeros(size(walkerHistM));
        ListDeadWalkersM =  zeros(size(walkerHistM));
        
        %% Distribute jobs (individuals)
        parfor iInd = 1:maximumThreadt
            distWalker = walkerHistM(:,iInd);
            NewWalkers = zeros(length(distWalker),1);
            DeadWalkers = zeros(length(distWalker),1);
            for iK = 1:length(distWalker)
                switch distWalker{iK}.status
                    case 1
                        [newWalker,dauWalker,deadWalker] = UpdatingIndivMobParCNBioCrustWating2(distWalker{iK},Vu,DeltaT,m,n,Lp);
                        distWalker{iK} = newWalker;
                        NewWalkers(iK) = dauWalker;
                        DeadWalkers(iK) = deadWalker;
                    case 2
                        distWalker{iK}.age = distWalker{iK}.age + DeltaT;
                        %distWalker{iK}.waitingT = distWalker{iK}.waitingT + DeltaT;
                end
            end
            walkerHistM(:,iInd) = distWalker(:);
            ListNewWalkersM(:,iInd) = NewWalkers(:);
            ListDeadWalkersM(:,iInd) = DeadWalkers(:);
        end
        
        %% Collect information of new borns and deads
        walkerHist = reshape(walkerHistM,1,numel(walkerHistM));
        tempD = reshape(ListDeadWalkersM,1,numel(walkerHistM));
        tempL = reshape(ListNewWalkersM,1,numel(walkerHistM));
        
        %% decayed cells return to the nutrients
        deads = find(tempD);
        if isempty(deads) == 0
            for i = 1:length(deads)
                index = deads(i);
                testWalker = walkerHist{index};
                position = testWalker.positionS;
                if testWalker.sp <5 %phototrophs
                    sitesCI(position(1),position(2),1) = sitesCI(position(1),position(2),1) + testWalker.V*rho/waterVolume(position(1),position(2)); %Dead body : cell lysis become the inactive cell and decay with a certain rate
                else %non-phototrophs
                    sitesCI(position(1),position(2),2) = sitesCI(position(1),position(2),2) + testWalker.V*rho/waterVolume(position(1),position(2)); %Dead body : cell lysis become the inactive cell and decay with a certain rate
                end
            end
        end
        
        for ii =1:6 %feedback to carbohydrates and ammonium / certain changes in biocarbonate
            for i = 1:2
                sitesC{ii}(:,:) = sitesC{ii}(:,:) +YDecay(i,ii)*sitesCI(:,:,i)*(1-exp(-decayRate(i)*DeltaT));
            end
        end
        
        for i = 1:2
            sitesCI(:,:,i) = sitesCI(:,:,i)*exp(-decayRate(i)*DeltaT);
        end
        
        %EPS degradation as a function of concentration of surrounding sugar? sitesC2{1}*exp(-decayRate(i)*DeltaT);
        % Hydrolysis of EPS
        HydrolysisEPS = sitesC2{1}.*(1-exp(-decayRate(3)*DeltaT*sitesC2{1}./EPSCcrit));
        sitesC2{1} = sitesC2{1} - HydrolysisEPS;
        sitesC{4} = sitesC{4} + HydrolysisEPS;
        
        
        %% New borns
        mothers = find(tempL);
        if isempty(mothers) == 0
            for i = 1:length(mothers)
                index = mothers(i);
                numberOfWalkers = numberOfWalkers + 1;
                testWalker = walkerHist{index};
                dauWalker1 = GenerateDaughterParCNBioCrust(testWalker,numberOfWalkers);
                walkerHist{numberOfWalkers} = dauWalker1;
                numberOfWalkers = numberOfWalkers + 1;
                dauWalker2 = GenerateDaughterParCNBioCrust(testWalker,numberOfWalkers);
                walkerHist{numberOfWalkers} = dauWalker2;
            end
        end
        
        %% for video :: for ensemble averages - you can remove it for single run
        PopulationMovie0 = zeros(m,n,numberOFP); % Selecting occupied patches for motility caculations
        for i = 1:numberOFP
            PopulationMovie0(:,:,i) = PopulationMovie{i,tt}(:,:);
            popWalkersPart(tt,i) = sum(sum(PopulationMovie{i,tt}(:,:)));
            dormpopWalkersPart(tt,i) = sum(sum(dormPopulationMovie{i,tt}(:,:)));
            muList{i}(:,:,tt) = muList{i}(:,:,tt)./popWalkersPart(tt,i);
        end
        
        
        
        %% Partition inorganic carbon following pH
        %sitesC = sitesCtemp;
        %sitesC2 = sitesC2temp;
        sitesC2t = sitesC{2};
        sitesC3 = sitesC{3};
        sitesC5 = sitesC{5};
        sitesC6 = sitesC{6};
        sitesC7 = sitesC{7};
        sitesC22 = sitesC2{2};
        sitesC23 = sitesC2{3};
        sitesC24 = sitesC2{4};
        sitesC25 = sitesC2{5};
        sitesC26 = sitesC2{6};
        sitesC27 = sitesC2{7};
        sitesCaTemp = sitesCa;
        parfor iInd = 1:m
            [C2t, C3, C5, C6, C7, C22, C23, C24, C25, C26, C27, CaTemp] = PHestimationMex(sitesC2t(iInd,:),sitesC3(iInd,:), sitesC5(iInd,:), sitesC6(iInd,:),sitesC7(iInd,:), sitesC22(iInd,:),sitesC23(iInd,:), sitesC24(iInd,:),sitesC25(iInd,:),sitesC26(iInd,:),sitesC27(iInd,:), sitesCaTemp(iInd,:,:), DeltaT, ZionM(iInd,:),pKTotList,thetaTC);
            sitesC2t(iInd,:) = C2t;
            sitesC3(iInd,:) = C3;
            sitesC5(iInd,:) = C5;
            sitesC6(iInd,:) = C6;
            sitesC7(iInd,:) = C7;
            sitesC22(iInd,:) = C22;
            sitesC23(iInd,:) = C23;
            sitesC24(iInd,:) = C24;
            sitesC25(iInd,:) = C25;
            sitesC26(iInd,:) = C26;
            sitesC27(iInd,:) = C27;
            sitesCaTemp(iInd,:,:) = CaTemp;
        end
        sitesC{2} = sitesC2t;
        sitesC{3} = sitesC3;
        sitesC{5} = sitesC5;
        sitesC{6} = sitesC6;
        sitesC{7} = sitesC7;
        sitesC2{2} = sitesC22;
        sitesC2{3} = sitesC23;
        sitesC2{4} = sitesC24;
        sitesC2{5} = sitesC25;
        sitesC2{6} = sitesC26;
        sitesC2{7} = sitesC27;
        sitesCa = sitesCaTemp;
        
        
        for i = 1:numberOFN
            sitesN{i,1} = sitesC{i,1}.*waterVolume;
            totalNutriPart{i,1}(tt) = sum(sum(sitesN{i,1}));
        end
        
    end
    
    %% Upload infomration of popualtion size and total nutrients with the time step of plottt
    for i = 1:numberOFP
        ActivityCells{i,examineT}(:,:) = muList{i}(:,:,plottt);
        Mumean(examineT,i) = mean(mean(muList{i}(:,:,plottt)));
        totalBioMass(examineT,i) = sum(sum(testV(:,:,i)))*rho;
        BioMassDist{i,examineT}(:,:) = testV(:,:,i)*rho;
        
        popWalkers(examineT,i) = popWalkersPart(plottt,i);
        dormpopWalkers(examineT,i) = dormpopWalkersPart(plottt,i);
        popMovie{i,examineT}(:,:) = PopulationMovie{i,plottt}(:,:);
        dormpopMovie{i,examineT}(:,:) = dormPopulationMovie{i,plottt}(:,:);
        measureMMovie{i,examineT}(:,:) = measVtot(:,:,i);
    end
    
    %timeConcDist2 = cell(5,1); %EPS, CO3, NH4, H+, PH
    for i = 1:numberOFN
        totNutConsumption{i,examineT}(:,:) = totCsum(:,:,i);
        GaseousEffluxMicrobes{i,examineT}(:,:) = totNgsum(:,:,i);
        GaseousEffluxBubbles{i,examineT}(:,:) = totNgBubble(:,:,i);
        totalNutri(examineT,i) = totalNutriPart{i,1}(plottt);
        timeConcDist{i,1}(:,:,examineT) = sitesC{i,1}(:,:);
    end
    
    timeConcDist{1,2}(:,:,examineT) = sitesCg{1,1}(:,:);
    timeConcDist{2,2}(:,:,examineT) = sitesCg{2,1}(:,:);
    timeConcDist{6,2}(:,:,examineT) = sitesCg{6,1}(:,:);
    timeConcDist{7,2}(:,:,examineT) = sitesCg{7,1}(:,:);
    
    for i = 1:7
        timeConcDist2{i,1}(:,:,examineT) = sitesC2{i};
    end
    
    timeConcDist2{8}(:,:,examineT) = sitesCa(:,:,1);
    timeConcDist2{9}(:,:,examineT) = sitesCa(:,:,2);
    
    sitesCaN1 = sitesCa(:,:,1).*waterVolume;
    sitesCaN2 = sitesCa(:,:,2).*waterVolume;
   
end

save(strcat('./',serial_id_new,'/',serial_idtemp,'/','Final.mat'),'-v7.3');
clear all;
toc;

