function [testWalker,dauWalkers,deadWalkers] = UpdatingIndivMobParCNBioCrustWating2(testWalker,Vu,dt1,m,n,Lp)
Vdmin = Vu*(2/1.433);
Vmin = Vdmin/5;
walkerID = testWalker.number;
muG = testWalker.muGcorr;
updateV = testWalker.V*(1+(muG*dt1)); %when delta t is defined as unit time
dauWalkers = 0;
deadWalkers = 0;
waitTime = testWalker.initT; %For inactivation
%reproduction
if muG<0
    
    %if waitTime > 0%Waiting Time 5 hours
    if waitTime > 24*60*60%Waiting Time 1 day
        
        if rand < 0.8 %probability of dying 20%
            testWalker.status = 2;
            testWalker.V = updateV;
            testWalker.age = testWalker.age + dt1;
            testWalker.waitingT = dt1;
            testWalker.initT = 0;
        else %Otherwise save the copy in the domain
            testWalker.status = 4;
            testWalker.V = updateV;
            deadWalkers = walkerID;
        end
        
    else if updateV < Vmin
            
            testWalker.status = 4;
            testWalker.V = updateV;
            deadWalkers = walkerID;
        else
            testWalker.V = updateV;
            testWalker.age = testWalker.age + dt1;
            testWalker.initT = testWalker.initT + dt1;
            %Calculate movement respected to growth rate of each species
            if testWalker.velocity ~= 0
                [testWalker] = MicrobeMobParCNBioCrust(testWalker,dt1,m,n,Lp);
            else
                jumpTProb = ones(7,1)/7;
                testWalker.fluxProb = testWalker.fluxProb + jumpTProb;
                testWalker.stayingT = testWalker.stayingT +dt1;
            end
            
        end
    end
    
else
    
    testWalker.initT = 0;
    if updateV >= Vdmin %big enough to generate daughters
        testWalker.status = 3;
        testWalker.V = updateV;
        dauWalkers = walkerID;
    else
        testWalker.V = updateV;
        testWalker.age = testWalker.age + dt1;
        %Calculate movement respected to growth rate of each species
        if testWalker.velocity ~= 0
            [testWalker] = MicrobeMobParCNBioCrust(testWalker,dt1,m,n,Lp);
        else
            jumpTProb = ones(7,1)/7;
            testWalker.fluxProb = testWalker.fluxProb + jumpTProb;
            testWalker.stayingT = testWalker.stayingT +dt1;
        end
    end
    
end

end
