potList = linspace(0.5,100,DryHlist);

for i = 1:length(potList)
    i
    pot1 = potList(i);
    pot = -pot1/9.81;
    parfor iInd = 1:m
        [WF,percol] = WaterDistAffineHTSatuAProfile(systemAffine(iInd,:,:),pot,Hurst,R1,R2,R);
        waterFilm(iInd,:) = WF(1,:);
        percolProb(iInd,:) = percol(1,:);
    end
    
    waterVolume = waterFilm*patchArea; %Volume  of water
    LocalWaterContentsList(:,:,i) = waterFilm./totaldM;
    LocalWaterSaturationList(:,:,i) = LocalWaterContents./porosityM;
    avgWC(i) = mean(mean(LocalWaterContentsList(:,:,i)));
    avgWS(i) = mean(mean(LocalWaterSaturationList(:,:,i)));
end