function [possibleDist, jumpTProbM] = MicrobeExpMobParBioCrustSpaceConstrain(HighDensityPatch,microbeVelocityM, normWaterM, gradApp, perShare)

global ed chi0 Lp

[m,n] = size(microbeVelocityM);
%      microbeVelocityM= (PopulationMovieT0>0).*microbeVelocityMT;
%      normWaterM = normaWaterT;
%      gradApp= gradApparent{i,1};
%      perShare = perShareT;
chi01 = chi0/2;
possibleDist = cell(m,n);
jumpTProbM = cell(m,n);

for i = 1:m
    for j = 1:n
        
        jumpTProbM{i,j} = zeros(7,1);
        
        microbeV = microbeVelocityM(i,j);
        netDistgC1 = gradApp(:,i,j);
        normWater = normWaterM{i,j};
        normnetDistgC1 = norm(netDistgC1);
        normnetDistgC = abs(normnetDistgC1);
        tortuosity = perShare(i,j);
                
        if microbeVelocityM(i,j) ~=0
            
            if normnetDistgC ~= 0
                
                tempGrowthC = chi01/microbeV;
                tempDirection = tempGrowthC*ed*netDistgC1;
                tempProb = exp(tempDirection).*normWater;
                test = isinf(tempProb);
                jumpTProb = zeros(7,1);
                if sum(test) == 0;
                    total = sum(tempProb);
                    if total ~= 0
                        jumpTProb = tempProb/sum(tempProb);
                    else
                        jumpTProb(7) = 1;
                    end
                else
                    [ctemp, cI] = max(abs(tempDirection).*normWater);
                    jumpTProb(cI) = 1;
                end
                
                possibleDisplacements = (microbeV*tortuosity.*jumpTProb)'*ed;
                
                possibleDist{i,j} = possibleDisplacements;
                jumpTProbM{i,j} = jumpTProb;
                
            else
                
                possibleDist{i,j} = [0 0];
                jumpTProbM{i,j} = normWater/sum(normWater);
                
            end
            
        else
            possibleDist{i,j} = [0 0];
            jumpTProbM{i,j} = zeros(7,1);
        end
        
        if rand<(HighDensityPatch(i,j)-1)
            possibleDist{i,j} = [1 0];
            temp = ones(7,1);
            temp(7) = 0;
            temp = temp.*normWater;
            jumpTProbM{i,j} = temp/sum(temp);
        end
        
    end
end
end

