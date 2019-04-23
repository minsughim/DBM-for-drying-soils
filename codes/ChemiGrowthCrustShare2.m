function [apparentGrowth5, apparentGrowth6, apparentGrowth7,apparentGrowth8] = ChemiGrowthCrustShare2(sitesC,sitesC2,MumaxT,nutShare)
global Ks mumax

%KpH = 9.04E-8;
KpH = 0.00001;
%KpH = 1;
KIhonoAOB = 0.5036*10; %Roughly estimated from the figure of Park 2009 modelling
KINH3AOB = 0.0729*10;  %Roughly estimated from the figure of Park 2009 modelling
KIhonoNOB = 0.5036; %Blackburnea 2007 
KINH3NOB = 0.0729;  %Blackburnea 2007

sitesH= sitesC2{4};

pHfeed = KpH./(KpH+sitesH);

compO2 = nutShare(:,:,1).*sitesC{1}./(Ks(5,1)+sitesC{1});
compSugar = nutShare(:,:,4).*sitesC{4}./(Ks(5,4)+sitesC{4});
compHCO3 = nutShare(:,:,3).*sitesC{3}./(Ks(5,3)+sitesC{3});
compNH3 = nutShare(:,:,6).*sitesC{6}./(Ks(5,6)+sitesC{6});
temp1 = min(compO2,compSugar);
temp1 = min(temp1,compNH3);
temp1 = min(temp1,compHCO3);
apparentGrowth5 = mumax(5).*pHfeed.*MumaxT.*temp1;

compO2 = (-1*Ks(6,1)./(-Ks(6,1)+sitesC{1}));
compSugar = nutShare(:,:,4).*sitesC{4}./(Ks(6,4)+sitesC{4});
compNO3 = nutShare(:,:,5).*sitesC{5}./(Ks(6,5)+sitesC{5});
temp1 = min(compO2,compSugar);
temp1 = min(temp1,compNO3);
apparentGrowth6 = mumax(6).*pHfeed.*MumaxT.*temp1;

compO2 = nutShare(:,:,1).*sitesC{1}./(Ks(7,1)+sitesC{1});
compCO2 = nutShare(:,:,2).*sitesC{2}./(Ks(7,2)+sitesC{2});
compHCO3 = nutShare(:,:,3).*sitesC{3}./(Ks(7,3)+sitesC{3});
compNH3 = nutShare(:,:,6).*sitesC{6}./(Ks(7,6)+sitesC{6});
compNH32 = KINH3AOB./(KINH3AOB+sitesC2{3});
conmpHNO2 = (KIhonoAOB./(KIhonoAOB+sitesC2{6}));
temp1 = min(compO2,compHCO3);
temp1 = min(temp1,compCO2);
temp1 = min(temp1,compNH3);
temp1 = min(temp1,compNH32);
temp1 = min(temp1,conmpHNO2);
apparentGrowth7 = mumax(7).*pHfeed.*MumaxT.*temp1;

compO2 = nutShare(:,:,1).*sitesC{1}./(Ks(8,1)+sitesC{1});
compHCO3 = nutShare(:,:,3).*sitesC{3}./(Ks(8,3)+sitesC{3});
compNO2 = nutShare(:,:,7).*sitesC{7}./(Ks(8,7)+sitesC{7});
compNH3 = KINH3NOB./(KINH3NOB+sitesC2{3});
conmpHNO2 = (KIhonoNOB./(KIhonoNOB+sitesC2{6}));
temp1 = min(compO2,compHCO3);
temp1 = min(temp1,compNO2);
temp1 = min(temp1,compNH3);
temp1 = min(temp1,conmpHNO2);
apparentGrowth8 = mumax(8).*pHfeed.*MumaxT.*temp1;

end
