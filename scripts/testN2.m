pause(30)
clear all
maxT = 20;
poolobj = parpool('local',maxT);

examineDay = 5;
examineDay2 = 5;
pot1 = 3;
plottt = 5;
maximumThreadt = maxT;

avgT = 25;
newT = 25;
desiccationIndex = 13;
lightOn2 = 0;
NH3ppb = 5;
HONOppb = 1;

for indexS = 5:8
tic
Main_BSC_photoY_drying3_WaterShortDecay_NH3_HONO_Temperature(NH3ppb, HONOppb, lightOn2,examineDay,desiccationIndex,avgT, newT, pot1, plottt, indexS,maximumThreadt)
toc
end

delete(poolobj)
