function [systemAffineT] = GenerateAffineSurface(meanPhi, H)
global m n dl

%distribution of percolation threshold (C.Du, C. Satik, and Y.C. Yotsos
%1996 : Percolation in a Fractional Browian Motion lattice)
mu1 = (H==1)*0.386 + (H==0.8)*0.39 + (H==0.65)*0.405 + (H==0.5)*0.42 + (H==0.35)*0.44 + (H==0.2)*0.46 + (H==0.1)*0.48 +(H==0)*0.5;
sigma1 = (H==1)*0.13 + (H==0.8)*0.12 + (H==0.65)*0.115+(H==0.5)*0.11+(H==0.35)*0.1+(H==0.2)*0.08+(H==0.1)*0.04 +(H==0)*0;
%pd = makedist('Normal','mu', mu1,'sigma',sigma1);

systemAffineT =zeros(m,n,2);

field1=Brownian_field(H,m,n,dl*m);
temp1 = min(min(field1));
translate1 = field1 - temp1+0.1;
a = meanPhi/mean(mean(translate1));
affinePhi = a*translate1;
D = 3-H;

for i= 1:m
    for j = 1:n
        systemAffineT(i,j,1) = affinePhi(i,j)*(affinePhi(i,j)<=1)+ 1*(affinePhi(i,j)>1); %fratal porosity
        systemAffineT(i,j,2) = D;
        affinePhitemp =0;
        while (affinePhitemp>0)*(affinePhitemp<1)==0
            affinePhitemp = normrnd(mu1, sigma1);
        end
        systemAffineT(i,j,3) = affinePhitemp; %p_c
        
    end
end

end
