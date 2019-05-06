pause(30)
maxT = 16;
poolobj = parpool('local',maxT);

examineDay = 5;
amountCations = 0.5;
avgT = 25;
reducedLevelC = 0.1;
pot1 = 3;
plottt = 5;
indexS = 1;
maximumThreadt = maxT;
for indexS = 1:3
    Main_BSC_noPG_noBio_noPH(examineDay,amountCations, avgT, reducedLevelC,pot1, plottt, indexS,maximumThreadt)
end
delete(poolobj)

%delete(gcp('nocreate'))
