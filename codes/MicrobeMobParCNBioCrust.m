function [Walker] = MicrobeMobParCNBioCrust(Walker,dt,m,n,Lp)

distance = norm(Walker.positionV);
tempPos= Walker.positionS;
y = tempPos(1);
x = tempPos(2);

if distance > Lp
    
    totalFluxProb = Walker.fluxProb/sum(Walker.fluxProb);
    newFluxProb = zeros(3);
    newFluxProb(2,2) = totalFluxProb(7);
    newFluxProb(2,3) = totalFluxProb(3);
    newFluxProb(2,1) = totalFluxProb(6);
    if rem(y,2) == 0
        newFluxProb(1,1) = totalFluxProb(1);
        newFluxProb(1,2) = totalFluxProb(2);
        newFluxProb(3,2) = totalFluxProb(4);
        newFluxProb(3,1) = totalFluxProb(5);
    else
        newFluxProb(1,2) = totalFluxProb(1);
        newFluxProb(1,3) = totalFluxProb(2);
        newFluxProb(3,3) = totalFluxProb(4);
        newFluxProb(3,2) = totalFluxProb(5);
    end
    
    MotilityB = reshape(newFluxProb',1,9);
    temp = randsample(1:9, 1, true, MotilityB);
    positiondY = ceil(temp(1)/3)-2;
    positiondX = rem(temp(1),3)-2;
    if positiondX == -2
        positiondX = 1;
    end
    %%%Boundary condition :: Periodic boundary conditions (x-direction.
    %%%Cylindirical)
    
    newX = mod(x+positiondX,n);
    newX = (newX==0)*n + (newX~=0)*newX; %periodic boundary conditions for column
    
    if y+positiondY==0
        newY = 1;
    elseif y+positiondY > m
        newY = m;
    else
        newY = y+positiondY;
    end
        
    if newX == x && newY == y
        newX = x;
        newY = y;
        Walker.stayingT = Walker.stayingT + dt;
    else
        Walker.stayingT = 0;
        Walker.velocity = 0;
        Walker.fluxProb = zeros(7,1);
        Walker.positionV = [0,0];
        Walker.mufield =zeros(2,1);
        Walker.normW = zeros(3);
    end
    
else
    newX = x;
    newY = y;
    Walker.stayingT = Walker.stayingT + dt;
end

Walker.positionS = [newY,newX];

end