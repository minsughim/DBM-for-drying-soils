pause(30)
clear all
maxT = 44;
poolobj = parpool('local',maxT);

examineDay = 5;
dryInd = 4;
pot1 = 3;
plottt = 5;
maximumThreadt = maxT;
pCO2frac = 1;
amountCations = 0.1;
avgT = 25;
alphaCO = 3;
lightOn =0;
NH3ppb = 50;
HONOppb = 0.005;
for indexS = 2:2
Main_BSC_photoY_drying3_H2ONO_100mu(NH3ppb, HONOppb, lightOn,examineDay,dryInd,  pot1, plottt, indexS,maxT)
end


delete(poolobj)
