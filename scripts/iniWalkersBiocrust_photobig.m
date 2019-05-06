%Random innoculation
%innoculSize = 5;
%posistionX1 = floor(0.5*m);
numberOfWalkers = 1;
PopulationMovie0 = zeros(m,n,numberOFP);
%Inoculation of autotrophs on the toplayer
%topM = floor(m/4);
for xx = 2:(m-1)
    for yy = 1:n
        popPH = initPopPH(xx,yy);
        for iW = 1:popPH
            for i = 1:4
                PopulationS(xx,yy,1) = PopulationS(xx,yy,1)+ 1;
                PopulationMapS{xx, yy,1}(PopulationS(xx,yy,1)) = numberOfWalkers;
                %Information of walkers
                iniWalker.number = numberOfWalkers;
                iniWalker.sp = i;
                iniWalker.V = Vu;
                iniWalker.age = 0;
                iniWalker.status = 1;
                iniWalker.positionS = [xx,yy];
                iniWalker.positionV = [0,0];
                iniWalker.velocity = microbeVelocityM(xx,yy);
                iniWalker.travelD = 0;
                iniWalker.stayingT = 0;
                iniWalker.waitingT = 0;
                iniWalker.initT = 0;
                iniWalker.fluxProb = zeros(7,1);
                iniWalker.mufield = zeros(2,1);
                iniWalker.normW = zeros(7,1);
                iniWalker.muGcorr = 0;                
                iniWalker.pg = pGmax;                
                walkerHist{numberOfWalkers} = iniWalker;
                
                numberOfWalkers = numberOfWalkers+1;
                
                PopulationMovie0(xx,yy,i) = PopulationMovie0(xx,yy,i) + 1;
                
            end
            
        end
    end
end

initPopph = sum(initPopPH(:));
inocPopothers = ceil(initPopph/(m-2)/n);
for xx = 2:(m-1)
    for yy = 1:n
        for iW = 1:inocPopothers
            for i = 5:numberOFP
                PopulationS(xx,yy,1) = PopulationS(xx,yy,1)+ 1;
                PopulationMapS{xx, yy,1}(PopulationS(xx,yy,1)) = numberOfWalkers;
                %Information of walkers
                iniWalker.number = numberOfWalkers;
                iniWalker.sp = i;
                iniWalker.V = Vu;
                iniWalker.age = 0;
                iniWalker.status = 1;
                iniWalker.positionS = [xx,yy];
                iniWalker.positionV = [0,0];
                iniWalker.velocity = microbeVelocityM(xx,yy);
                iniWalker.travelD = 0;
                iniWalker.stayingT = 0;
                iniWalker.waitingT = 0;
                iniWalker.initT = 0;
                iniWalker.fluxProb = zeros(7,1);
                iniWalker.mufield = zeros(2,1);
                iniWalker.normW = zeros(7,1);
                iniWalker.muGcorr = 0;
                iniWalker.pg = pGmax; 
                walkerHist{numberOfWalkers} = iniWalker;                
                numberOfWalkers = numberOfWalkers+1;
                PopulationMovie0(xx,yy,i) = PopulationMovie0(xx,yy,i) + 1;
                
            end
            
        end
    end
end


numberOfWalkers = numberOfWalkers - 1;