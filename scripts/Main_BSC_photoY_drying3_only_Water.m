function Main_BSC_photoY_drying3_only_Water(NH3ppb, HONOppb, lightOn2,examineDay,desiccationIndex, avgT, newT, pot1, plottt, indexS,maximumThreadt)
rng('shuffle')
dryingT = 1; % 24 hours of drying for all in this code 

global m n Rad3 Rad3Over2 dt Lp dl patchArea  R1 R2 %modelling parameters
global V0 chi0 mumax Ks Ymax mrate Vu R rho % prameters for microbes
global dataPoints ed listE numberOFN numberOFP testMinT Hurst
global averageT amplitudeT attenuationRate MaximumIncidence repsiRatio
global iniConcentrationGas iniConcentration

pCO2frac = 1;
amountCations = 0.1;
% avgT = 25;
peneD = 0.2;
% pot1 = 3;
% plottt = 5;
% indexS = 9;
alphaCO = 3;

%serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);
serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,pot1,pCO2frac,amountCations,avgT,alphaCO, indexS);

cd(serial_id_new);
serial_idtemp_dd = sprintf('Dinural_cycle');
cd(serial_idtemp_dd);
load('Final.mat')
cd ..
cd ..

%% Drying
examineDay = dryingT;
reducedLevelC = 10^(-16); %to start with zero input from bottom, it increases slowly depending on the hydration condition
satuIndex = 0; % saturation index 1: indicating completely saturated case, saturation index 0: indicating field capacit corresponds to 0.66 saturation condition
repsiRatio = 0.1; %
Zratio = amountCations; %nut pH = 7.5;
Zion = Zratio*100*10^(-3)/40; %in mol/L (net ion concentration that t
calciumF = amountCations;
initCalcium = calciumF*100*10^(-3)/40; %%100 maicro g/g

%% Time information for the fixed time step
dt = 60; %second %1 min.
testMinT = 0.1; % second: diffusion equation resolution of 1 second.

examineH = dryingT*24 + 24*5; %sustain as the dry soil one more day
timeLine = (dt*plottt)/3600:(dt*plottt)/3600:examineH;
TinitalD = 4*(3600/dt)*24; %examine time, (seconds)*hours
T = (3600/dt)*examineH; %examine time, (seconds)*hours
Td = (3600/dt)*dryingT*24; %examine time, (seconds)*hours
dataPoints = T/plottt;
timeLapseList = zeros(dataPoints,1);
DeltaT = dt;
preDryHlist = TinitalD/plottt;
dryPot = 50;
DryHlist = Td/plottt;
load('Dessication_Biocrust.mat', 'potListTotal')
potList(1:DryHlist) = potListTotal(:,desiccationIndex);

%% Assinging the property of soil necessary
%water film thickness and matric potentialfor different depth.
pot1 = dryPot;
pot = -pot1/9.81;
[waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,Hurst,R1,R2,R);
[microbeVelocityM, probOccup] = velocityMicrobeMatrix2(waterFilm, pot,V0,R);
perShare = ones(m,n);
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

ListUnsatu = find(waterSatu<1);
if isempty(ListUnsatu)==1
    exRatewithoutDiff =zeros(m,n);
else
    exRatewithoutDiff = ((localDiffG./(1-LocalWaterSaturation))-(Aqdiff./LocalWaterSaturation)).*sepecifInterA./TotalporeD;
end

tic

%% Starting dynamics
timeM = zeros(m,n);

waterContentsDynamics = zeros(m,n,dataPoints);
waterVolumeDynamics = zeros(m,n,dataPoints);
waterFilmDynamics = zeros(m,n,dataPoints);
InvasedIslandDynamics = zeros(m,n,dataPoints);
gasVolumeDynamics = zeros(m,n,dataPoints);

for examineT = 1:dataPoints
    
    examineT
    
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
    end
    
    
    avgWF(examineT) = mean(waterFilm(:));
    avgWC(examineT) = mean(LocalWaterContents(:));
    avgWS(examineT) = mean(LocalWaterSaturation(:));
    avgGC(examineT) = mean(LocalGasContents(:));
    
    waterContentsDynamics(:,:,examineT) = LocalWaterContents;
    waterVolumeDynamics(:,:,examineT) = waterVolume;
    waterFilmDynamics(:,:,examineT) = waterFilm;
    InvasedIslandDynamics(:,:,examineT) = InvasedIsland;
    gasVolumeDynamics(:,:,examineT) = gasVolume;
    
end

serial_id3 = sprintf('WaterDynamics_Drying%d_varyingT',desiccationIndex);
save(strcat('./',serial_id_new,'/',serial_id3,'.mat'),'waterContentsDynamics', 'waterVolumeDynamics', 'waterFilmDynamics', 'InvasedIslandDynamics','gasVolumeDynamics');

clear all;
toc;

