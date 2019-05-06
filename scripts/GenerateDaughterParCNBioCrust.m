function dauWalker = GenerateDaughterParCNBioCrust(testWalker,numberOfWalkers);
 

dauWalker = testWalker;
dauWalker.number = numberOfWalkers;
dauWalker.V = (testWalker.V)/2;
dauWalker.age = 0;
dauWalker.status = 1;
dauWalker.velocity = testWalker.velocity;
dauWalker.positionS = testWalker.positionS;
dauWalker.positionV = testWalker.positionV;
dauWalker.travelD = 0;
dauWalker.stayingT = 0;
dauWalker.waitingT = 0;
dauWalker.initT = 0;
dauWalker.fluxProb = zeros(7,1);
dauWalker.mufield = zeros(2,1);
dauWalker.normW = zeros(3);
dauWalker.muGcorr = 0;

end