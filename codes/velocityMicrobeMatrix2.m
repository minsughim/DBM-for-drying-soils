function [microbeVelocityM, percolProb] = velocityMicrobeMatrix2(filmthickness, pot,V0,R);
%global V0 R
[m,n] = size(filmthickness);
rca = 7.41*10^(-6)/(-pot);
microbeVelocityM = zeros(m,n);

for i = 1:m
    for j = 1:n
        
        d = filmthickness(i,j);
        
        ratio = R/d;
        lambdap = 1/(1- 9/16*ratio +1/8*ratio^3);
        lambdan = 1+ 9/8*ratio +(9/8*ratio)^2;
        lambda = sqrt(lambdap^2 + lambdan^2);
        M_force=6*3.14*0.001*R*V0;
        
        if d>2*R
            microbeVelocityM(i,j) = V0/lambda;
        else
            
            dataH=rca*(2*R-2*d)/(rca+R);
            d=(2*d+dataH)/2;
            AA=rca+2*d-R;
            BB=R+rca;
            cos_theta=(2*AA*BB+(4*AA^2*BB^2-4*BB^2*AA^2)^0.5)/BB^2/2;
            sin_theta=(1-cos_theta^2)^0.5;
            force=(2*pi*R*sin_theta*0.072*sin_theta-pi*(R*sin_theta)^2*pot)*1e-2;
            
            if force<0
                force=0;
            end;
            microbeVelocityM(i,j) = V0*((1/lambda)-(force/M_force));%restricted velocity\
            
            if microbeVelocityM(i,j) <=0
                microbeVelocityM(i,j) =0;
            end;
            %       else
            %           microbeVelocityM(i,j) = 0;
            %       end
        end
    end
end

veloMatric = microbeVelocityM;

results2 = zeros(2,1);

[realIsland,town2,results2,numberOfCluster] = IslandStatHex(veloMatric,0);

temp2 = size(results2);

if  temp2(1) == 1
    
    largestIsland = results2(2);
    
else if length(results2) == 0
        
        largestIsland = 0;
        
    else
        largestIsland = max(results2(:,2));
        
        
    end
end

percolProb = largestIsland/(m*n);



end
