function [systemAffine] = GenerateAffineSurfaceHM(meanPhi, H)
global m n

%distribution of percolation threshold (C.Du, C. Satik, and Y.C. Yotsos
%1996 : Percolation in a Fractional Browian Motion lattice)
mu1 = (H==1)*0.386 + (H==0.8)*0.39 + (H==0.65)*0.405 + (H==0.5)*0.42 + (H==0.35)*0.44 + (H==0.2)*0.46 + (H==0.1)*0.48 +(H==0)*0.5;
sigma1 = (H==1)*0.13 + (H==0.8)*0.12 + (H==0.65)*0.115+(H==0.5)*0.11+(H==0.35)*0.1+(H==0.2)*0.08+(H==0.1)*0.04 +(H==0)*0;
%pd = makedist('Normal','mu', mu1,'sigma',sigma1);

systemAffine =zeros(m,n,2);
D = 3-H;
systemAffine(:,:,1) = meanPhi*ones(m,n);
systemAffine(:,:,2) = D*ones(m,n);
systemAffine(:,:,3) = mu1*ones(m,n); %p_c

end
