function [apparentGrowth,respiGrowth] = PhotoGrowthShare(sitesC, sitesC2,MumaxTt,intensityProfile,nutShare)
global Ks mumax repsiRatio

KsP = 1.4;
KsI = 295;
%KpH = 9.04E-8;
KpH = 0.00001;
%KpH = 1;
[m,n] = size(sitesC{1});
clear photoGrowth respiGrowth apparentGrowth
respiGrowth(m,n,8) =0; 
apparentGrowth(m,n,4) =0;

sitesH = sitesC2{4};

pHfeed = KpH./(KpH+sitesH);

muIntensity = intensityProfile./(intensityProfile + KsP + intensityProfile.*intensityProfile/KsI);

tempRespi = min(nutShare(:,:,1).*sitesC{1}./(Ks(1,1)+sitesC{1}),nutShare(:,:,6).*sitesC{6}./(Ks(2,6)+sitesC{6})); %O2 and TAN
tempRespi = min(tempRespi,nutShare(:,:,4).*sitesC{4}./(Ks(1,4)+sitesC{4})); %sugar 
tempRespi = min(tempRespi,nutShare(:,:,3).*sitesC{3}./(Ks(3,3)+sitesC{3})); %HCO3- 
tempRespi = min(tempRespi,KsP./(KsP+intensityProfile)); %inhibit by light intensity
respiGrowth1 = repsiRatio*mumax(1).*MumaxTt.*tempRespi;

for i = 1:4 %all phototrohps

    compCO2 = (Ks(i,2)>=0).*nutShare(:,:,2).*(sitesC{2}./(Ks(i,2)+sitesC{2})) + (Ks(i,2)<0).*(-1*Ks(i,2)./(-Ks(i,2)+sitesC{2}));
    compHCO3 = (Ks(i,3)>=0).*nutShare(:,:,3).*(sitesC{3}./(Ks(i,3)+sitesC{3})) + (Ks(i,3)<0).*(-1*Ks(i,3)./(-Ks(i,3)+sitesC{3}));
    compNO3 = (Ks(i,5)>=0).*nutShare(:,:,5).*(sitesC{5}./(Ks(i,5)+sitesC{5})) + (Ks(i,5)<0).*(-1*Ks(i,5)./(-Ks(i,5)+sitesC{5}));
    compNH3 = (Ks(i,6)>=0).*nutShare(:,:,6).*(sitesC{6}./(Ks(i,6)+sitesC{6})) + (Ks(i,6)<0).*(-1*Ks(i,6)./(-Ks(i,6)+sitesC{6}));
    %[compCO2(1,1), compHCO3(1,1), compNO3(1,1), compNH3(1,1)]    
    temp1 = min(compCO2,compHCO3);
    temp1 = min(temp1,compNO3);
    temp1 = min(temp1,compNH3);
    apparentGrowth(:,:,i) = mumax(1).*pHfeed.*MumaxTt.*min(temp1,muIntensity);    
    respiGrowth(:,:,i) = respiGrowth1.*pHfeed;       
end


end