%parallel.importProfile('/cluster/apps/matlab/support/EulerLSF8h.settings')
%cluster = parcluster('EulerLSF8h')
%cluster.SubmitArguments = '-W 72:00 -R "rusage[mem=5000]"'
maxT = 8;
poolobj = parpool('local',maxT);

examineDay =5;
examineDay2 = 3;
NfixationR = 0.1;
amountCations = 0.5;
avgT = 25;
peneD = 0.2;
reducedLevelC = 0.1;
pot1 = 3;
plottt = 5;
indexS = 1;
maximumThreadt = maxT;
tic
Main_BSC_noPG_wetting(examineDay,examineDay2,NfixationR, amountCations, avgT, peneD,reducedLevelC,pot1, plottt, indexS,maximumThreadt)
toc
delete(poolobj)

%delete(gcp('nocreate'))
