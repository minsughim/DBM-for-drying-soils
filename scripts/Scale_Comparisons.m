rng('shuffle')

global m n Rad3 Rad3Over2 dt Lp dl patchArea  R1 R2 %modelling parameters
global dataPoints ed listE numberOFN numberOFP testMinT Hurst

desiccationIndex = 6;

load('Dessication_Biocrust.mat', 'potListTotal')
potList = potListTotal(:,desiccationIndex);

m = 1;
n = 1;
R1 = 10^(-4);
R2 = 10^(-7);
R = 1e-6; % radius of the bacteria : assume that its shape is cylindrical
Rad3Over2 = sqrt(3)/2;
Rad3 = sqrt(3);
Lp = 0.00005/Rad3Over2;
dl = Lp/Rad3; % length of the edge of hexagon (top layer)
patchArea = Rad3Over2*dl*dl;

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

D0 = 10^(-4)/(24*3600);

%% Assinging the property of soil necessary
%water film thickness and matric potentialfor different depth.

WFlist = zeros(length(potList),1);
GFlist = zeros(length(potList),1);
WSlist = zeros(length(potList),1);
WClist = zeros(length(potList),1);
GClist = zeros(length(potList),1);
AqDifflist = zeros(length(potList),1);
GasDifflist = zeros(length(potList),1);
specAlist = zeros(length(potList),1);
exRatelist = zeros(length(potList),1);


tic
for i = 1:length(potList)
    i
%for i = 1:5  
    pot1 = potList(i);
    pot = -pot1/9.81;
    
    [waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,Hurst,R1,R2,R);
    waterVolume = waterFilm*patchArea; %Volume  of water
    TotalporeD = capilaryPored.*1.5;
    totaldM = TotalporeD./porosity;
    LocalWaterContents = waterFilm./totaldM;
    LocalWaterSaturation = LocalWaterContents./porosity;
    LocalGasContents = porosity-LocalWaterContents;
    Aqdiff = D0*(LocalWaterContents.^(10/3)).*porosity.^(-2);
    localDiffG = D0*10^4*(LocalGasContents.^(10/3)).*porosity.^(-2);
    gasVolume = LocalGasContents.*totaldM*patchArea; %Volume of water
    
    if waterSatu == 1
        exRatewithoutDiff =zeros(m,n);
    else
        exRatewithoutDiff = ((localDiffG./(1-LocalWaterSaturation))-(Aqdiff./LocalWaterSaturation)).*sepecifInterA./TotalporeD;
    end

    WFlist(i) = waterFilm;
    GFlist(i) = LocalGasContents.*totaldM;
    WSlist(i) = LocalWaterSaturation;
    WClist(i) = LocalWaterContents;
    GClist(i) = LocalGasContents;
    AqDifflist(i) = Aqdiff;
    GasDifflist(i) = localDiffG;
    specAlist(i) = sepecifInterA;
    exRatelist(i) = exRatewithoutDiff;
       
end
toc
