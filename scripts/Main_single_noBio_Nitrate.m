function [sitesC2ini, sitesCaini, sitesCgini, sitesCini, waterFilm, LocalGasContents] = Main_single_noBio_Nitrate(PartialPressure,amountCations,nitrateF,avgT, pot1)
rng('shuffle')
%To calcauted initial condition for the steady state condition
examineDay = 50;
%pot1 = 3;
plottt = 5;
initMoles = 0.01;
dryingT = 0;
frac = 0.9;

%Zratio = 0.5; %nut pH = 7.5;
Zratio = amountCations; %nut pH = 7.5;
Zion = Zratio*100*10^(-3)/40; %in mol/L (net ion concentration that t
calciumF = amountCations;
initCalcium = calciumF*100*10^(-3)/40; %%100 micro g/g in mole

%lambda = 8.8*1e-6; %Thermal conductivity
averageT = avgT; %[degree C]
amplitudeT = 0; %[degree C]
attenuationRate = 0.4*0.001; %[mm^(-1)]
MaximumIncidence = 500; %[?mol.m?2.s?1]

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
%%%%%Target
%m = 300; %Number of rows for top layer
%n = 50; %Number of columns for top layer
%Lp = 0.0001; %[m]
%%%%%
m = 1; %Number of rows for top layer
n = 1; %Number of columns for top layer
LayerDepth = 0.000;
depthList = linspace(0,LayerDepth,m);
%Lp = (depthList(2)-depthList(1))/Rad3Over2; % size of patch for toplayer [m]
Lp = 0.0002/Rad3Over2;
dl = Lp/Rad3; % length of the edge of hexagon (top layer)
R1 = 10^(-4);
R2 = 10^(-7);
R = 1e-6; % radius of the bacteria : assume that its shape is cylindrical

patchArea = Rad3Over2*dl*dl;

%% Time information for the fixed time step
%dt = 60; %second %1 min.
dt = 120; %second %1 min.
testMinT = 0.1; % second: diffusion equation resolution of 1 second.
examineH = dryingT + 24*examineDay;
timeLine = (dt*plottt)/3600:(dt*plottt)/3600:examineH;
T = (3600/dt)*examineH; %examine time, (seconds)*hours
Td = (3600/dt)*dryingT; %examine time, (seconds)*hours
dataPoints = T/plottt;
timeLapseList = zeros(dataPoints,1);
DeltaT = dt;
DryHlist = Td/plottt;
potseries = 3*ones(dataPoints,1);
%potseries(end-DryHlist+1:end) = logspace(log10(pot1),1,DryHlist);
potseries(end-DryHlist+1:end) = linspace(pot1,100,DryHlist);


%% Parameters (biological and chemical)
numberOFN = 8; %O2, CO2, HCO3-, sugar, NO3-, NH4+, NO2- N2O: used for food or byproduct
numberOFP = 8; % four phototrphs(nitrogen fixing, subindexed), one aerobe, one anaerobe (denitrifying), one chemoautotrophic AOB, and NOB.
KA = 10^(-9.4003);
K2C = 10^(-6.3819);
K1C = 10^(-10.3767);

pa = 1.01325e5;
density_air = air_density(averageT, 0,pa)*1000;

%PartialPressure= [pO2,pCO2,pNH3,pHONO,pN2O,pN2];
%% Atmospheric input O2, CO2, N2
iniConcentrationGasO2 = density_air*PartialPressure(1); %mg/L
iniConcentrationGasCO2 = density_air*PartialPressure(2); %mg/L #383ppm
iniConcentrationGasN2 = density_air*PartialPressure(9); %mg/L
iniConcentrationGas = zeros(numberOFN,1);
iniConcentrationGas(1) = iniConcentrationGasO2; %mg/L
iniConcentrationGas(2) = iniConcentrationGasCO2; %mg/L
iniConcentrationGas(6) = density_air*PartialPressure(6); %mg/L  a few ppb: assumed 5 %Gong et al 2011: atmospheric level of ammonia
iniConcentrationGas(7) = density_air*PartialPressure(7); %mg/L %1ppb for HONO
iniConcentrationGas(8) = density_air*5*PartialPressure(8); %mg/L %500ppb (o.5ppm) for N2O
% Dry composition information http://www.engineeringtoolbox.com/air-composition-d_212.html
%% All dissolved %chemical vaues from johnson2005relevance
iniConcentration(1) = iniConcentrationGasO2*0.0318; % N4 : O2 mg/L : Henry's law at 1 atm
iniConcentration(2) = iniConcentrationGasCO2*0.8553; % N1 : CO2 mg/L : Henry's law at 1 atm
iniConcentration(3) = iniConcentration(2)*61.0171/44.0098*K1C/10^(-7.5); % at equilbrium with atmospheric level
iniConcentration(4) = 26.9958*initMoles; % N5 : SS 1M (readily degradable organic substrate) : CH1.5O0.5N0.1 %Initmole only controls the sugar concentraion. (Catillo-Monroy 2010): 34.61mg/kg soil: with porosity 0.3: 34.61*2.6*0.3~
iniConcentration(5) = nitrateF; % N2 : NO3- : H (Abed 2013, for hot desert with lots of nitrate)
iniConcentration(6) = 0.88*18/14; %johnson2005relevance
iniConcentration(7) = 2*(1-frac); %NO2- same as NO3- to start with
iniConcentration(8) = iniConcentrationGas(8)*0.6110; %N2O

%% Assinging the property of soil necessary
%water film thickness and matric potentialfor different depth.
%pot1 = 3;
pot = -pot1/9.81;
meanPhi = 0.6;
H = 0.35;
porosity = 0.4; %Field capacity at water saturation 1/2
%Homogeneous profile%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1 = (H==1)*0.386 + (H==0.8)*0.39 + (H==0.65)*0.405 + (H==0.5)*0.42 + (H==0.35)*0.44 + (H==0.2)*0.46 + (H==0.1)*0.48 +(H==0)*0.5;
sigma1 = (H==1)*0.13 + (H==0.8)*0.12 + (H==0.65)*0.115+(H==0.5)*0.11+(H==0.35)*0.1+(H==0.2)*0.08+(H==0.1)*0.04 +(H==0)*0;
systemAffine =zeros(m,n,2);
D = 3-H;
systemAffine(:,:,1) = meanPhi*ones(m,n);
systemAffine(:,:,2) = D*ones(m,n);
systemAffine(:,:,3) = mu1*ones(m,n); %p_c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Heterogeneous profile%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,H,R1,R2,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

porosityM = porosity*ones(m,n);
perShare = 1;
waterVolume = waterFilm*patchArea; %Volume  of water
TotalporeD = capilaryPored.*1.5;
totaldM = TotalporeD./porosityM;
LocalWaterContents = waterFilm./totaldM;
LocalWaterSaturation = LocalWaterContents./porosityM;
LocalGasContents = porosityM-LocalWaterContents;
Aqdiff = (LocalWaterContents.^(10/3)).*porosityM.^(-2);
localDiffG = 10^4*(LocalGasContents.^(10/3)).*porosityM.^(-2);
gasVolume = LocalGasContents.*totaldM*patchArea; %Volume of water
InvasedIsland = ones(m,n);
ListUnsatu = find(LocalWaterSaturation<1);
if isempty(ListUnsatu)==1
    exRatewithoutDiff =zeros(m,n);
else
    exRatewithoutDiff = ((localDiffG./(1-LocalWaterSaturation))-(Aqdiff./LocalWaterSaturation)).*sepecifInterA./TotalporeD;
end

%% Temperature conditions and Light intensity
lambdaMatrix = ((1-porosityM)/2.9 + LocalWaterContents/0.57 + LocalGasContents/0.0012).^(-1);%Soil Thermal conductivity [W/mK]
CvMatrix = 10^(6)*(1.94*(1-porosityM)+4.189*LocalWaterContents+ 0.0012*LocalGasContents); %Soil Thermal diffusivity[m^2/s]
alphaMatrix = lambdaMatrix./CvMatrix; %Soil Thermal diffusivity[m^2/s]

%% Preallocate variables

totalNutriPart = cell(numberOFN,1);
tempsitesC = cell(numberOFN,1);
changeNurient = cell(numberOFN,1);

sitesC = cell(numberOFN,1); %concentration of nutrients
sitesN = cell(numberOFN,1); % mass of nutrients
sitesCg = cell(numberOFN,1); %concentration of nutrients
sitesNg = cell(numberOFN,1); %concentration of nutrients
sitesC2 = cell(7,1); %concentration of auxiliary nutrients %EPS, CO3, NH3, H+, PH HONO Ca
sitesN2 = cell(7,1); %concentration of auxiliary nutrients %EPS, CO3, NH3, H+, PH HONO Ca


reactionC = cell(numberOFN,1);
timeConcDist = cell(numberOFN,1);
timeConcDist2 = cell(5,1); %EPS, CO3, NH3, H+, PH HONO Ca
totalNutri = zeros(dataPoints,numberOFN);

%%N2O or N2 production can be included later
HenryConstList = cell(numberOFN,1);
HenryConstProfile = cell(numberOFN,1);
NT = cell(numberOFN,1);


%% calcauting environmental conditions and preallocate variables
%Only z-depth vary : needs to be changed when the physical domain is
%heterogeneous
%No henry's law applied for nitrogen. (denitrification is not included
%yet : to consider HONO or NO2, you have to include henry's law for nitrogen gases)
[intensityList,temperatureList,MumaxT, HenryConstList{1,1},HenryConstList{2,1},HenryConstList{6,1},HenryConstList{7,1},HenryConstList{8,1}, pKList] = EnvironmentProfile2(depthList,timeLine,alphaMatrix,averageT, amplitudeT, attenuationRate, MaximumIncidence);
for i = 3:5
    HenryConstList{i,1} = zeros(dataPoints,m);
end

% modify the concentraiton of nitrate based on the amount of nitrate
nitrateNamount = TotalporeD*patchArea*nitrateF;
iniConcentration(5) = nitrateNamount/waterVolume;

for i = 1:numberOFN%except sugar and EPS (initial condition for only CO2 and O2)
    sitesC{i,1} = iniConcentration(i)*ones(m,n); % Initial condition of concentrations are based on the aqueous phase
    sitesCg{i,1} = iniConcentrationGas(i);
    changeNurient{i,1} = zeros(m+2, n+2);
    totalNutriPart{i,1} = zeros(plottt,1);
    timeConcDist{i,1} = zeros(m, n,dataPoints);
    timeConcDist{i,2} = zeros(m, n,dataPoints);
    if i < 8
        timeConcDist2{i,1} = zeros(m, n,dataPoints);
    end
end


sitesC2{1} = zeros(m,n); %EPS
sitesC2{2} = 0; %CO3
sitesC2{3} = iniConcentrationGas(6)*1490.8*ones(m,n); %NH3 (almost zero)
sitesC2{4} = 10^(-7.5)*ones(m,n); %[H+] in M
sitesC2{5} = 7.5*ones(m,n); %pH
sitesC2{6} = iniConcentrationGas(7)*1221.9*ones(m,n); %HONO
sitesC2{7} = initCalcium*ones(m,n); %Ca2+ (calcium ion) in M (Johnson 2005 relevance)
% Calcite
sitesCa = zeros(m,n,2); % calcite complex / precipitation
sitesCa(:,:,1) = sitesC2{2}.*sitesC2{7}*10^(3.22);
sitesCa(:,:,2) = sitesC2{2}.*sitesC2{7}/10^(-8.48);

for i = 1:8
    sitesN{i} = sitesC{i}.*waterVolume;
end

for i = 1:7
    sitesN2{i} = sitesC2{i}.*waterVolume;
end

effluxList = zeros(dataPoints,5);
ListGas = [1 2 6 7 8];
AlphaM = cell(numberOFN,1);
sitesNEquil = cell(numberOFN,1);
HenryDomain = cell(numberOFN,1);
temp2 = zeros(m,n);
for iG = 1:length(ListGas)
    gI = ListGas(iG);
    for i = 1:n
        temp2(:,i) = HenryConstList{gI,1}(1,:);
    end
    HenryDomain{gI} = temp2;
    AlphaM{gI} = temp2.*waterVolume + gasVolume;
    sitesNEquil{gI} = temp2.*iniConcentrationGas(gI);
end


waterSatuT = zeros(dataPoints,1);
WCT = zeros(dataPoints,1);
waterVolumeT = zeros(dataPoints,1);
gasVolumeT = zeros(dataPoints,1);


for examineT = 1:dataPoints
    
    if examineT > dataPoints-DryHlist +1
        
        pot1 = potseries(examineT);
        pot = -pot1/9.81;
        
        [waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,Hurst,R1,R2,R);
        waterVolume = waterFilm*patchArea; %Volume  of water
        LocalWaterContents = waterFilm./totaldM;
        LocalWaterSaturation = LocalWaterContents./porosityM;
        LocalGasContents = porosityM-LocalWaterContents;
        Aqdiff = (LocalWaterContents.^(10/3)).*porosityM.^(-2);
        localDiffG = 10^4*(LocalGasContents.^(10/3)).*porosityM.^(-2);
        gasVolume = LocalGasContents.*totaldM*patchArea; %Volume of water
        
        %% For boundary conditions
        InvasedIsland = 1;
        exRatewithoutDiff = ((localDiffG./(1-LocalWaterSaturation))-(Aqdiff./LocalWaterSaturation)).*sepecifInterA./TotalporeD;
        lambdaMatrix = ((1-porosityM)/2.9 + LocalWaterContents/0.57 + LocalGasContents/0.0012).^(-1);%Soil Thermal conductivity [W/mK]
        CvMatrix = 10^(6)*(1.94*(1-porosityM)+4.189*LocalWaterContents+ 0.0012*LocalGasContents); %Soil Thermal diffusivity[m^2/s]
        alphaMatrix = lambdaMatrix./CvMatrix; %Soil Thermal diffusivity[m^2/s]
        
    end
    waterSatuT(examineT) = LocalWaterSaturation;
    WCT(examineT) = LocalWaterContents;
    waterVolumeT(examineT) = waterVolume;
    gasVolumeT(examineT) = gasVolume;
    
    [intensityList,temperatureList,MumaxT, HenryConstList{1,1},HenryConstList{2,1},HenryConstList{6,1},HenryConstList{7,1},HenryConstList{8,1}, pKList] = EnvironmentProfile2(depthList,timeLine,alphaMatrix,averageT, amplitudeT, attenuationRate, MaximumIncidence);
    
    for i = 1:numberOFN%except sugar and EPS (initial condition for only CO2 and O2)
        if i < 8
            sitesC2{i} = sitesN2{i}./waterVolume;
        end
        sitesC{i} = sitesN{i}./waterVolume;
    end
    
    for iG = 1:length(ListGas)
        gI = ListGas(iG);
        for i = 1:n
            temp2(:,i) = HenryConstList{gI,1}(1,:);
        end
        HenryDomain{gI} = temp2;
        AlphaM{gI} = temp2.*waterVolume + gasVolume;
        sitesNEquil{gI} = temp2.*iniConcentrationGas(gI);
    end
    %% water informaion update
    for i = 1:numberOFN
        totalNutriPart{i,1} = zeros(plottt,1);
    end

    changeN = zeros(m,n,numberOFN);
    
    for tt = 1:plottt
        
        t = (examineT-1)*plottt + tt;
    
        for i = 1:numberOFN
            sitesN{i,1} = waterVolume.*sitesC{i,1};
            if i < 8
                sitesN2{i,1} = waterVolume.*sitesC2{i,1};
            end
        end
        for i = 1:numberOFN+3
            reactionC{i,1} = zeros(m,n);
        end
        
        
        %% Reset Environmental conditions (Light and temperature -> growth rate accordingly)
        thetaTC(:,:) = diag(temperatureList(1,:))*ones(m,n);
        pKTotList = zeros(m,n,6); %pK
        for i = 1:6
            pKTotList(:,:,i) = diag(pKList(1,:,i))*ones(m,n);
        end
        
        %calculate the efflux and equilibrisind the gas-liquid phase: mass
        %conserve
        [efflux, sitesC, sitesCg, sitesC2] = EquilibriumConcentrationFlux2(iniConcentrationGas,waterVolume,gasVolume,sitesC,sitesC2,sitesCg,HenryDomain,InvasedIsland,AlphaM,sitesNEquil,changeN);
        effluxList(examineT,:) = effluxList(examineT,:) + efflux;% + sum(sum((timeConcDist{1,2}(:,:,examineT)-sitesCg{1,1}.*gasVolume.*InvasedIsland)));
        
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
        ZionM  = Zion;
        
        
        for iInd = 1:n
           [C2t, C3, C5, C6, C7, C22, C23, C24, C25, C26, C27, CaTemp] = PHestimationMex(sitesC2t(:,iInd),sitesC3(:,iInd), sitesC5(:,iInd), sitesC6(:,iInd),sitesC7(:,iInd), sitesC22(:,iInd),sitesC23(:,iInd), sitesC24(:,iInd),sitesC25(:,iInd),sitesC26(:,iInd),sitesC27(:,iInd), sitesCaTemp(:,iInd,:), DeltaT, ZionM(:,iInd),pKTotList,thetaTC);
            
            sitesC2t(:,iInd) = C2t;
            sitesC3(:,iInd) = C3;
            sitesC5(:,iInd) = C5;
            sitesC6(:,iInd) = C6;
            sitesC7(:,iInd) = C7;
            sitesC22(:,iInd) = C22;
            sitesC23(:,iInd) = C23;
            sitesC24(:,iInd) = C24;
            sitesC25(:,iInd) = C25;
            sitesC26(:,iInd) = C26;
            sitesC27(:,iInd) = C27;
            sitesCaTemp(:,iInd,:) = CaTemp;
            
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
            if i <8
                sitesN2{i,1} = sitesC2{i,1}.*waterVolume;
            end
        end
        
    end
    
    %% Upload infomration of total nutrients with the time step of plottt
    
    %timeConcDist2 = cell(5,1); %EPS, CO3, NH4, H+, PH
    for i = 1:numberOFN
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

sitesC2ini = sitesC2;
sitesCaini = sitesCa;
sitesCgini = sitesCg;
sitesCini = sitesC;

%figure(1)
%temp = zeros(dataPoints,1);
%temp(:) = timeConcDist2{5};
%hold on
%plot(timeLine,temp)

%figure(2)
%hold on
%totVerticalArea = mean(totaldM(:))*Lp*n/porosity;
%plot(timeLine,effluxList*10^6/dt/totVerticalArea,'DisplayName','effluxList')
%figure(2)
%semilogy(timeLine,effluxList(:,4)*3600/dt/plottt/totVerticalArea,'DisplayName','effluxList')
%hold on
%semilogy(timeLine,-1*effluxList(:,4)*3600/dt/plottt/totVerticalArea,'DisplayName','effluxList')

%save(strcat('./',serial_id,'/','Final.mat'),'-v7.3');
%clear all;
%toc;

