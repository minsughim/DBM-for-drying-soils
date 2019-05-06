pause(60)
maxT = 16;
poolobj = parpool('local',maxT);

examineDay = 5;
amountCations = 0.1;
avgT = 25;
reducedLevelC = 0.1;
pot1 = 3;
plottt = 5;
indexS = 1;
maximumThreadt = maxT;
for indexS = 1:3
    Main_BSC_par_noBio(examineDay,amountCations, avgT, reducedLevelC,pot1, plottt, indexS,maximumThreadt)
end
delete(poolobj)
