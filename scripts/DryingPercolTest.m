
global m n Rad3 Rad3Over2 dt Lp dl patchArea  R1 R2 %modelling parameters
global V0 chi0 mumax Ks Ymax mrate Vu R rho % prameters for microbes
global dataPoints ed listE numberOFN numberOFP testMinT
global averageT amplitudeT attenuationRate MaximumIncidence repsiRatio
global iniConcentrationGas iniConcentration


examineDay = 1;
pot1 = 3;
plottt = 5;
amountCations = 0.1;
avgT = 25;
pCO2frac = 1;
indexS = 0;
initMoles = 1;
maximumThreadt = 2;

reducedLevelC = 0.1;
NfixationR = 0.1;
peneD = 0.2;
repsiRatio = 0.1; %

Zratio = amountCations; %nut pH = 7.5;
Zion = Zratio*100*10^(-3)/40; %in mol/L (net ion concentration that t
calciumF = amountCations;
initCalcium = calciumF*100*10^(-3)/40; %%100 maicro g/g
satuIndex = 0; % saturation index 1: indicating completely saturated case, saturation index 0: indicating field capacit corresponds to 0.66 saturation condition
germT = 1; % 1 hour: response delay rather than the germination time

if satuIndex == 0 % At field condition (water Saturation 0.66)-> not completely saturated case (amount of water held by capilarity only)
    serial_id = sprintf('HT_BSC_Day%dPot%.1fpCO2frac%dCat%.1favT%dindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,indexS);
else
    serial_id = sprintf('HT_BSC_Satu_Day%dPot%.1fpCO2frac%dCat%.1favT%dindex%d',examineDay,pot1, pCO2frac,amountCations,avgT,indexS);
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
%%%%%Target
%m = 300; %Number of rows for top layer
%n = 50; %Number of columns for top layer
%Lp = 0.0001; %[m]
%%%%%
m = 300; %Number of rows for top layer
n = 100; %Number of columns for top layer
LayerDepth = 0.02;
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
timeLapseList = zeros(dataPoints,1);
DeltaT = dt;

%% Parameters (biological and chemical)

%partial pressure of trace gases

pN2 = 0.7809;
pO2 = 0.2095;
pCO2 = 0.000383*pCO2frac; %383 ppm
pNH3 = 5*10^(-9); %5ppb: Gong et al 2011 %2ppb in McCally
pHONO = 10^(-9); %1ppb for HONO : Su et al 2011 from dry area
pN2O = 5*10^(-7); %mg/L %500ppb (0.5ppm) for N2O
PartialPressure =[pO2,pCO2,0,0,0,pNH3,pHONO,pN2O,pN2];
[sitesC2ini, sitesCaini, sitesCgini, sitesCini] = Main_single_noBio(PartialPressure,amountCations,avgT);

Parameters_CNratio_BSC_netlogo


%% Temperature conditions and Light intensity
%lambda = 8.8*1e-6; %Thermal conductivity
averageT = avgT; %[degree C]
amplitudeT = 5; %[degree C]
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
[systemAffine] = GenerateAffineSurface(meanPhi, Hurst);
%[systemAffine] = GenerateRandomSurface(meanPhi, Hurst);
%[waterFilm,percolProb,waterSatu,Aqdiff,totalPored] = WaterDistAffineHTSatu(systemAffine, pot,Hurst);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
porosityM = (2/3)*porosity*(ones(m,n)+ rand(m,n));
potWet = 0.5;
pot = -potWet/9.81;
[waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,Hurst,R1,R2,R);
TotalporeD = capilaryPored.*(1+rand(m,n));
totaldM = TotalporeD./porosityM;

reducedLevelC = 0;
potList = linspace(potWet,5,10);
ListRealIsland = zeros(m,n,50);
ListInvasedIsland = zeros(m,n,50);
ListWaterContents = zeros(m,n,50);
ListwaterSaturation = zeros(m,n,50);

for iP = 1:50
    iP
    pot1 = potList(iP);
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
    %Aqdiff = ones(m,n);
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
    ListRealIsland(:,:,iP) = realIsland;
    ListInvasedIsland(:,:,iP) = InvasedIsland;
    ListWaterContents(:,:,iP) = LocalWaterContents;
    ListwaterSaturation(:,:,iP) = LocalWaterSaturation;
    
    
end