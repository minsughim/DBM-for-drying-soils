%pause(60)
clear all
maxT = 44;
poolobj = parpool('local',maxT);

examineDay = 1;
pot1 = 3;
plottt = 5;
NH3ppb = 5;
HONOppb = 1;
indexS = 1;

for indexS = 1:1
Main_BSC_noBio_H2ONO(examineDay, pot1, NH3ppb, HONOppb,plottt, indexS,maxT)
end


delete(poolobj)
