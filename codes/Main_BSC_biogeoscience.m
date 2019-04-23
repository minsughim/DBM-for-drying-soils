function Main_BSC_biogeoscience(examineDay, pot1, plottt, indexS)
rng('shuffle')

global m n Rad3 Rad3Over2 dt Lp dl patchArea  R1 R2 %modelling parameters
global V0 chi0 mumax Ks Ymax mrate Vu R rho % prameters for microbes
global dataPoints ed listE numberOFN numberOFP testMinT
global averageT amplitudeT attenuationRate MaximumIncidence repsiRatio
global iniConcentrationGas iniConcentration


% These are the example inputs for the main code

%examineDay = 5; %Time for the dynamics (given in days)
%pot1 = 3; % matric potential for unsaturated soils, with the unit of [-kPa]
%plottt = 5; % This value is for saving. dt*plottt is the time interval. As dt is given as 60 sec below, plottt = 5 would save variables for every 5 mins.
%indexS = 0; %index for ensemble averages
maximumThreadt = feature('numcores'); %number of cores for the parallel computing

if pot1 == 0 % for the completly saturated case
    serial_id = sprintf('HT_BSC_satu_Day%d_index%d',examineDay, indexS);
    satuIndex = 1;
else
    serial_id = sprintf('HT_BSC_Day%d_Pot%.1f_index%d',examineDay,pot1, indexS); 
    satuIndex = 0;
end

mkdir(serial_id);

%% Spatial information

Rad3Over2 = sqrt(3)/2;
Rad3 = sqrt(3);
ed(1,:) = [-0.5,Rad3Over2];
ed(2,:) = [0.5, Rad3Over2];
ed(3,:) = [1,0];
ed(4,:) = [0.5,-1*Rad3Over2];
ed(5,:) = [-0.5, -1*Rad3Over2];
ed(6,:) = [-1,0];
ed(7,:) = [0,0];

m = 300; %Number of rows for top layer
n = 20; %Number of columns for top layer
LayerDepth = 0.02; % 2cm
depthList = linspace(0,LayerDepth,m);
Lp = (depthList(2)-depthList(1))/Rad3Over2; % size of patch for toplayer [m]
dl = Lp/Rad3; % length of the edge of hexagon (top layer)
R1 = 10^(-4);
R2 = 10^(-7);
patchArea = Rad3Over2*dl*dl;

g=graph('2d-deg6-trilattice',{n,m,Lp});
listV = g.vertices;
listE = g.edges;
%periodic boundary condition for y direction
listE(m*n-1,m) = 1;
listE(m*n,m) = 1;
listE(m,m*n-1) = 1;
listE(m,m*n) = 1;
for i = 1:2:(m-1)
    listE((n-1)*m+i ,i) = 1;
    listE(i,(n-1)*m+i) = 1;
end
for i = 2:2:(m-2)
    listE((n-1)*m+i-1 ,i) = 1;
    listE((n-1)*m+i ,i) = 1;
    listE((n-1)*m+i+1 ,i) = 1;
    listE(i,(n-1)*m+i-1) = 1;
    listE(i,(n-1)*m+i) = 1;
    listE(i,(n-1)*m+i+1) = 1;
end

%% Time information for the fixed time step
dt = 60; %second %1 min.
testMinT = 0.1; % second: diffusion equation resolution of 1 second.
examineH = 24*examineDay;
%timeLine = 0:(dt*plottt)/3600:examineH;
timeLine = (dt*plottt)/3600:(dt*plottt)/3600:examineH;
T = (3600/dt)*examineH; %examine time, (seconds)*hours
dataPoints = T/plottt;
DeltaT = dt;

%% Parameters (biological and chemical)

% The amount of cations relates the soil pH and the set-up point for the
% local pH values 
amountCations = 0.1; % amount of background cation in soils (uniformly distributed)
Zratio = amountCations; %nut pH = 7.5 in case 1;
Zion = Zratio*100*10^(-3)/40; %in mol/L (net ion concentration that t
calciumF = amountCations; 
initCalcium = calciumF*100*10^(-3)/40; %%100 micro g/g

% Different scenarios for the environmental conditions can be considered
% here by altering, CO2 fraction or average temperature
pCO2frac = 1; % for the scenario of higher CO2 fraction
avgT = 25; % the diurnal cycles of temperature is given with a sin function. Here, the avearge is given as 25 degree celcius.
reducedLevelC = 0.1; % boundary condition for the nutrients

%partial pressure of trace gases
pN2 = 0.7809;
pO2 = 0.2095;
pCO2 = 0.000383*pCO2frac; %383 ppm
pNH3 = 5*10^(-9); %5ppb: Gong et al 2011 %2ppb in McCally
pHONO = 10^(-9); %1ppb for HONO : Su et al 2011 from dry area
pN2O = 5*10^(-7); %mg/L %500ppb (0.5ppm) for N2O
PartialPressure =[pO2,pCO2,0,0,0,pNH3,pHONO,pN2O,pN2];
[sitesC2ini, sitesCaini, sitesCgini, sitesCini] = Main_single_noBio(PartialPressure,amountCations,avgT);

Parameters_CNratio_BSC_netlogo_photoY2

%% Temperature conditions and Light intensity
%lambda = 8.8*1e-6; %Thermal conductivity
averageT = avgT; %[degree C]
amplitudeT = 5; %[degree C]
peneD = 0.2; % penetration depth for the light
attenuationRate = peneD*0.001; %[mm^(-1)]
MaximumIncidence = 500; %[?mol.m?2.s?1]

%% Assinging the property of soil necessary
%water film thickness and matric potentialfor different depth.
%pot1 = 3;
pot = -pot1/9.81;
meanPhi = 0.6;
Hurst = 0.35;
porosity = 0.4; %Field capacity at water saturation 1/2
%Homogeneous profile%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%systemAffine = GenerateAffineSurfaceHM(meanPhi, Hurst);
%[waterFilm,percolProb,waterSatu,Aqdiff,totalPored,sepecifInterA] = WaterDistAffineHMSatuA(systemAffine, pot,Hurst);
%[microbeVelocityM, probOccup] = velocityMicrobeMatrixHM(waterFilm, pot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Heterogeneous profile%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[systemAffine] = GenerateAffineSurface(meanPhi, Hurst);
[systemAffine] = GenerateRandomSurface(meanPhi, Hurst);
%[waterFilm,percolProb,waterSatu,Aqdiff,totalPored] = WaterDistAffineHTSatu(systemAffine, pot,Hurst);
[waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,Hurst,R1,R2,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

porosityM = (2/3)*porosity*(ones(m,n)+ rand(m,n));

if satuIndex == 0 % At field condition (water Saturation 0.66)-> not completely saturated case (amount of water held by capilarity only)
    [microbeVelocityM, probOccup] = velocityMicrobeMatrix2(waterFilm, pot,V0,R);
    perShare = percolProb*probOccup;
    waterVolume = waterFilm*patchArea; %Volume  of water
    TotalporeD = capilaryPored.*(1+rand(m,n));
    totaldM = TotalporeD./porosityM;
    LocalWaterContents = waterFilm./totaldM;
    LocalWaterSaturation = LocalWaterContents./porosityM;
    LocalGasContents = porosityM-LocalWaterContents;
    Aqdiff = (LocalWaterContents.^(10/3)).*porosityM.^(-2);
    localDiffG = 10^4*(LocalGasContents.^(10/3)).*porosityM.^(-2);
    gasVolume = LocalGasContents.*totaldM*patchArea; %Volume of water
    
    %% For boundary conditions
    [realIsland,town2,results2,numberOfCluster] = IslandStatHex(LocalGasContents,0.2); %Invasive percolation for BCC case Pc = 0.2
    
    listInvasionBottom =realIsland(end,:);
    ListLocation = find(listInvasionBottom);
    InvasedIsland = zeros(m,n);
    for i = 1:length(ListLocation)
        IslandIndex = listInvasionBottom(ListLocation(i));
        InvasedIsland(realIsland==IslandIndex) = reducedLevelC;
    end
    InvasedIsland(end,:) = reducedLevelC;
    
    listInvasion =realIsland(1,:);
    ListLocation = find(listInvasion);
    for i = 1:length(ListLocation)
        IslandIndex = listInvasion(ListLocation(i));
        InvasedIsland(realIsland==IslandIndex) = 1;
    end
    InvasedIsland(1,:) = 1;
    
    ListUnsatu = find(waterSatu<1);
    if isempty(ListUnsatu)==1
        exRatewithoutDiff =zeros(m,n);
    else
        exRatewithoutDiff = ((localDiffG./(1-LocalWaterSaturation))-(Aqdiff./LocalWaterSaturation)).*sepecifInterA./TotalporeD;
    end
    
else %satuIndex == 1 completly saturated case (biocrust soaked in water)
    [microbeVelocityM, probOccup] = velocityMicrobeMatrix2(waterFilm, pot,V0,R);
    perShare = ones(m,n);
    TotalporeD = capilaryPored.*(1+1.2*rand(m,n));
    totaldM = TotalporeD./porosityM;
    waterVolume = TotalporeD*patchArea; %Volume  of water
    LocalWaterContents = porosityM;
    LocalWaterSaturation = ones(m,n);
    LocalGasContents = zeros(m,n);
    Aqdiff = (LocalWaterContents.^(10/3)).*porosityM.^(-2);
    localDiffG = zeros(m,n);
    gasVolume = zeros(m,n); %Volume of water
    
    %% For boundary conditions
    InvasedIsland = zeros(m,n);
    InvasedIsland(1,:) = 1;
    exRatewithoutDiff =zeros(m,n);
end

MaxCapacityPop = capilaryPored.*patchArea/Vdmin;

lambdaMatrix = ((1-porosityM)/2.9 + LocalWaterContents/0.57 + LocalGasContents/0.0012).^(-1);%Soil Thermal conductivity [W/mK]
CvMatrix = 10^(6)*(1.94*(1-porosityM)+4.189*LocalWaterContents+ 0.0012*LocalGasContents); %Soil Thermal diffusivity[m^2/s]
alphaMatrix = lambdaMatrix./CvMatrix; %Soil Thermal diffusivity[m^2/s]

save(strcat('./',serial_id,'/','Parameters.mat'));

%% Preallocate variables
PopulationS = zeros(m,n);
dormPopulationS = zeros(m,n);
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
sitesC = cell(numberOFN,1); %concentration of nutrients
sitesN = cell(numberOFN,1); % mass of nutrients
sitesCg = cell(numberOFN,1); %concentration of nutrients
sitesNg = cell(numberOFN,1); %concentration of nutrients
sitesC2 = cell(7,1); %concentration of auxiliary nutrients %EPS, CO3, NH3, H+, PH HONO Ca
ZionM = Zion*ones(m,n);

timeConcDist = cell(numberOFN,1);
timeConcDist2 = cell(5,1); %EPS, CO3, NH3, H+, PH HONO Ca
totalNutri = zeros(dataPoints,numberOFN);
diffMatrix = zeros(m,n,numberOFN);
diffMatrix2 = zeros(m,n,5);

%%N2O or N2 production is not included in this code
HenryConstList = cell(numberOFN,1);
HenryConstProfile = cell(numberOFN,1);
NT = cell(numberOFN,1);

producedSugar = zeros(m,n);
totCsum = zeros(m,n,numberOFN);
totNgsum = zeros(m,n,numberOFN);
changeN = zeros(m,n,numberOFN);
reactionC = zeros(m,n,numberOFN+3);
reactionCEPS = zeros(m,n);
PopulationS = zeros(m,n);
dormPopulationS = zeros(m,n);
testV = zeros(m,n,numberOFP);
testVinactive = zeros(m,n,numberOFP);
measVtot = zeros(m,n,numberOFP);
pop = zeros(m,n,numberOFP);
wakeUp = zeros(m,n,numberOFP);
muSp = zeros(m,n,numberOFP);
extraEPS = zeros(m,n);

ACtivityDepth = log(MaximumIncidence/sqrt(KsP*KsI))*attenuationRate;

initPopPH = zeros(m,n); % Initial distribution of phototrophs are assigned based on the photoactivity for faster equilibration
for i = 1:n
    initPopPH(:,i) = 30*exp(-(depthList-ACtivityDepth)/ACtivityDepth); %EPS
end
initPopPH = ceil(initPopPH);

%% Information for bubble formation
pN2endOthers = 1 - sum(PartialPressure(1:8)); %non-reactive gaseous elements
totNgBubble = zeros(m,n,numberOFN);
P_a = 1.4;

%% Inoculation with individual walker information
iniWalkersBiocrust_photo

%% calcauting environmental conditions and preallocate variables
%function of depth only (lateral difference of temperature is not included
%in this code)
[intensityList,temperatureList,MumaxT, HenryConstList{1,1},HenryConstList{2,1},HenryConstList{6,1},HenryConstList{7,1},HenryConstList{8,1}, pKList,Density_air] = EnvironmentProfileDensityAirTmax(depthList,timeLine,alphaMatrix);

for i = 3:5
    HenryConstList{i,1} = zeros(dataPoints,m);
end

rhoAir = zeros(m,n);
for i = 1:n
    rhoAir(:,i) = Density_air(1,:); %EPS
end

for i = 1:numberOFN %except sugar and EPS (initial condition for only CO2 and O2)
    sitesC{i,1} = iniConcentration(i)*ones(m,n); % Initial condition of concentrations are based on the aqueous phase
    sitesCg{i,1} = PartialPressure(i)*rhoAir;
    changeNurient{i,1} = zeros(m+2, n+2);
    totalNutriPart{i,1} = zeros(plottt,1);
    timeConcDist{i,1} = zeros(m, n,dataPoints);
    timeConcDist{i,2} = zeros(m, n,dataPoints);
    if i < 8
        sitesC2{i,1} = sitesC2ini{i}*ones(m,n);
        timeConcDist2{i,1} = zeros(m, n,dataPoints);
    end
end

sitesCa = zeros(m,n,2); %
sitesCa(:,:,1) = sitesCaini(1)*ones(m,n);
sitesCa(:,:,2) = sitesCaini(2)*ones(m,n);

sitesC2{1} = zeros(m,n); %Disrtirbution of EPS is pre-assgined following the light distribution (allows the system to equilibrate faster)
for i = 1:n
    sitesC2{1}(:,i) = EPSCcrit*exp(-(depthList-ACtivityDepth)/ACtivityDepth); %EPS
end

for i = 1:8
    sitesN{i} = sitesC{i}.*waterVolume;
end

for i = 1:7
    sitesN2{i} = sitesC2{i}.*waterVolume;
end

sitesCI = zeros(m,n,2);
decayRate = [mrate(1),mrate(5), 0.34/24/3600]; %EPS degradation rate from Wolf 2007

effluxList = zeros(dataPoints,5);
ListGas = [1 2 6 7 8]; %Gas that can bubble out (ebullition)
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

for examineT = 1:dataPoints
    
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
    
    for iG = 1:length(ListGas)
        gI = ListGas(iG);
        for i = 1:n
            temp2(:,i) = HenryConstList{gI,1}(examineT,:);
            rhoAir(:,i) = Density_air(examineT,:); %EPS
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
    
    totCsum = 0*totCsum;
    totNgsum = 0*totNgsum;
    changeN = 0*changeN ;
    reactionC = 0*reactionC;
    
    for tt = 1:plottt 
        
        t = (examineT-1)*plottt + tt;
        
        if rem(t, (12*60*60/dt)) == 0; %every 12 hours save essential information of microbial populations and substrate concentrations
            exposedHours = t*dt/(60*60)
            serial_id2 = sprintf('HTBioCrustCNcycleHour%d',exposedHours);
            save(strcat('./',serial_id,'/',serial_id2,'.mat'),'popWalkers', 'dormpopWalkers', 'popMovie','dormpopMovie','totalBioMass','timeConcDist','timeConcDist2','effluxList','timeLine','sitesC','sitesC2', 'sitesCa','totNutConsumption', 'GaseousEffluxMicrobes', 'Mumean','dataPoints', 'totalNutri','exposedHours','Lp','m','n','numberOFP','numberOFN','averageT','depthList','dt','plottt');
        end
        
        for i = 1:numberOFN
            sitesN{i,1} = waterVolume.*sitesC{i,1};
            if i < 8
                sitesN2{i,1} = waterVolume.*sitesC2{i,1};
            end
        end
        
        %% Reset Environmental conditions (Light and temperature -> calculating growth rates accordingly)
        intensityProfile = diag(intensityList(examineT,:))*ones(m,n);
        MumaxTt(:,:) = diag(MumaxT(examineT,:))*ones(m,n);
        thetaTC(:,:) = diag(temperatureList(examineT,:))*ones(m,n);
        
        pKTotList = zeros(m,n,6); %pK
        for i = 1:6
            pKTotList(:,:,i) = diag(pKList(examineT,:,i))*ones(m,n);
        end
        
        %calculate the efflux and equilibrisind the gas-liquid phase: mass
        %conservation
        [efflux, sitesC, sitesCg, sitesC2] = EquilibriumConcentrationDensity(waterVolume,gasVolume,sitesC,sitesC2,sitesCg,HenryDomain,InvasedIsland,AlphaM,sitesNEquil,sitesNGEquil,changeN);
        effluxList(examineT,:) = effluxList(examineT,:) + efflux;
        
        gradApparent = cell(numberOFP,2);
        possibleDist = cell(numberOFP,2);
        jumpTProbM = cell(numberOFP,2);
        wakeUp = 0*wakeUp;
        
        %% calcualting mumax based on the new environmental conditions
        [apparentGrowth(:,:,1:4),respiGrowthR] = PhotoGrowth(sitesC, sitesC2, MumaxTt,intensityProfile);
        [apparentGrowth(:,:,5),apparentGrowth(:,:,6),apparentGrowth(:,:,7),apparentGrowth(:,:,8)] = ChemiGrowthCrust2(sitesC,sitesC2,MumaxTt);
        
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
        
        %% Water might be not enough for the photosynthesis: Then shutdown photosynthesis -> for very dry case
        
       % reactionCwater = abs(reactionC(:,:,10)).*(reactionC(:,:,10)<0);
       % WaterShortage = (reactionCwater<waterVolume*10^(3)).*(1-reactionCwater./waterVolume*10^(-3)); %change in water volume larger than 0.1% than, consider it as water shortage.
        
        % reevaluate photogrowth term based on water shortage
       % for i = 1:4
       %     apparentGrowth(:,:,i) = WaterShortage.*apparentGrowth(:,:,i);
       % end
        
        % Nitrite oxidiser also requires water
       % apparentGrowth(:,:,8) = WaterShortage.*apparentGrowth(:,:,8);

        
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
        
        %% Diffusion Process: Diffusion-reaction process determines time step.        
        %EPS concentration changes the diffusivity of substrates in aqueous
        %phase
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
                    Diffresult{i} = DiffusionProcess_aq_gas(diffMatrix(:,:,i),waterFilm,waterVolume,gasVolume,sitesC{i}(:,:),sitesCg{i}(:,:),HenryDomain{i}(:,:),reactionC(:,:,i),InvasedIsland,sitesNEquil{i}(:,:),sitesNGEquil{i}(:,:),AlphaM{i}(:,:),listE,Lp,dt,testMinT, reducedLevelC);
                case 2 %CO2
                    Diffresult{i} = DiffusionProcess_aq_gas(diffMatrix(:,:,i),waterFilm,waterVolume,gasVolume,sitesC{i}(:,:),sitesCg{i}(:,:),HenryDomain{i}(:,:),reactionC(:,:,i),InvasedIsland,sitesNEquil{i}(:,:),sitesNGEquil{i}(:,:),AlphaM{i}(:,:),listE,Lp,dt,testMinT, reducedLevelC);
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
                    Diffresult{i} = DiffusionProcess_aq_gas(diffMatrix(:,:,i),waterFilm,waterVolume,gasVolume,sitesC{i}(:,:),sitesCg{i}(:,:),HenryDomain{i}(:,:),reactionC(:,:,i),InvasedIsland,sitesNEquil{i}(:,:),sitesNGEquil{i}(:,:),AlphaM{i}(:,:),listE,Lp,dt,testMinT, reducedLevelC);
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
                    Diffresult2{i} = DiffusionProcess_aq_gas(diffMatrix2(:,:,i),waterFilm,waterVolume,gasVolume,sitesC2{i}(:,:),sitesCg{6}(:,:),HenryDomain{6}(:,:),zeros(m,n),InvasedIsland,sitesNEquil{6}(:,:),sitesNGEquil{6}(:,:),AlphaM{6}(:,:),listE,Lp,dt,testMinT, reducedLevelC);
                case 4 %HONO
                    Diffresult2{i} = DiffusionProcess_aq_gas(diffMatrix2(:,:,i),waterFilm,waterVolume,gasVolume,sitesC2{6}(:,:),sitesCg{7}(:,:),HenryDomain{7}(:,:),zeros(m,n),InvasedIsland,sitesNEquil{7}(:,:),sitesNGEquil{7}(:,:),AlphaM{7}(:,:),listE,Lp,dt,testMinT, reducedLevelC);
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
        
        %% Ebullition when the built up partial pressure of gasous elements are larger than the fraction P_a
        
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
            limitF1 = min(nutShare(:))
            limitF2 = max(nutShare(:))
            DeltaT
            %Effect of pH applies to the next time step: Here only mass
            %share is applied
            [apparentGrowth(:,:,1:4),respiGrowthR] = PhotoGrowthShare(sitesC, sitesC2,MumaxTt,intensityProfile,nutShare);
            [apparentGrowth(:,:,5),apparentGrowth(:,:,6),apparentGrowth(:,:,7),apparentGrowth(:,:,8)] = ChemiGrowthCrustShare2(sitesC, sitesC2, MumaxTt,nutShare);
        end
       
        % reevaluate photogrowth term based on water shortage
       % for i = 1:4
       %     apparentGrowth(:,:,i) = WaterShortage.*apparentGrowth(:,:,i);
       % end
        
        % Nitrite oxidiser also requires water
       % apparentGrowth(:,:,8) = WaterShortage.*apparentGrowth(:,:,8);
       
        totalGrowth = apparentGrowth + respiGrowthR;
              
        %% Update possible displacement and jumping probability based on the grwoth rate field (microbial mibility)
        for i = 1:numberOFP
            gradApparent{i,1} = CalculateGradientField(totalGrowth(:,:,i),mumax(i),m,n,Lp);
            [possibleDist{i,1}, jumpTProbM{i,1}] = MicrobeExpMobParBioCrustSpaceConstrain(HighDensityPatch,(PopulationMovie0(:,:,i)>0).*microbeVelocityM, normWater, gradApparent{i,1}, perShare);
        end
        
        %% Update individuals based on spatial information
        
        muSp = 0*muSp;
        
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
        
        %% Distribute jobs (individuals): Parallel computing enters here:
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
        
        
        %% Newborns 
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
        
        %% for video :: for ensemble averages - you can remove it for a single run
        PopulationMovie0 = zeros(m,n,numberOFP); % Selecting occupied patches for motility caculations
        for i = 1:numberOFP
            PopulationMovie0(:,:,i) = PopulationMovie{i,tt}(:,:);
            popWalkersPart(tt,i) = sum(sum(PopulationMovie{i,tt}(:,:)));
            dormpopWalkersPart(tt,i) = sum(sum(dormPopulationMovie{i,tt}(:,:)));
            muList{i}(:,:,tt) = muList{i}(:,:,tt)./popWalkersPart(tt,i);
        end
        
        %% Partition inorganic carbon following pH
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
        %tic
        parfor iInd = 1:m %The pH calculation is the most time consuming part in this model after diffusion processes (Mex files are provided for Mac Os. and Linux only) 
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
        %toc
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
    
    
end

save(strcat('./',serial_id,'/','Final.mat'),'-v7.3');
clear all;
toc;

