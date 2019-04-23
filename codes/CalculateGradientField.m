function [gradApparent] = CalculateGradientField(apparentGrowthT,mumax,m,n,Lp);
global ed

gradApparent = zeros(2,m,n);
apparentGrowthTemp = generateBCcylinderY(apparentGrowthT);

for i = 1:m
    for j = 1:n
        
        apparentGrowthT = apparentGrowthTemp(i:i+2,j:j+2);
        
        gC13 = (apparentGrowthT(2,3)-apparentGrowthT(2,1))/Lp;
        
        if rem(i,2) == 0
            gC11 = (apparentGrowthT(1,1)-apparentGrowthT(3,2))/Lp;
            gC12 = (apparentGrowthT(1,2)-apparentGrowthT(3,1))/Lp;
        else
            gC11 = (apparentGrowthT(1,2)-apparentGrowthT(3,3))/Lp;
            gC12 = (apparentGrowthT(1,3)-apparentGrowthT(3,2))/Lp;
        end
        
        distgC1 = zeros(3,2);
        distgC1(1,:) = gC11*ed(1,:);
        distgC1(2,:) = gC12*ed(2,:);
        distgC1(3,:) = gC13*ed(3,:);
        
        netDistgC1 = sum(distgC1);
        gradApparent(:,i,j) = netDistgC1/mumax;
        
    end
end

end

