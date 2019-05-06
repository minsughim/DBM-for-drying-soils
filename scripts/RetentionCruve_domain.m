potList = linspace(3,50,100);
tic
for i = 1:1
    i   
    pot1 = -potList(i)/9.81;
    [waterFilm,percolProb,waterSatu,capilaryPored,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot1,Hurst,R1,R2,R);
    LocalWaterContents = waterFilm./totaldM;
    LocalWaterSaturation = LocalWaterContents./porosityM;
    LocalGasContents = porosityM-LocalWaterContents;
    Aqdiff = (LocalWaterContents.^(10/3)).*porosityM.^(-2);
    localDiffG = 10^4*(LocalGasContents.^(10/3)).*porosityM.^(-2);
   
    WaterContents(:,:,i) = LocalWaterContents;
    GasContents(:,:,i) = LocalGasContents;
    WaterSaturation(:,:,i) = LocalWaterSaturation;
    DiffusionAq(:,:,i) = Aqdiff;
    DiffusionGas(:,:,i) = localDiffG; 
end    
toc