pCO2frac = 1;
NH3ppb =  20;
HONOppb = 0.005;
amountCations = 0.1;
avgT = 25;
nitrateF = 5000;

pot1 = 1000;

pN2 = 0.7809;
pO2 = 0.2095;
pCO2 = 0.000383*pCO2frac; %383 ppm
pNH3 = NH3ppb*10^(-9); %default: 5ppb: Gong et al 2011 %2ppb in McCally
pHONO = HONOppb*10^(-9); %1ppb for HONO : Su et al 2011 from dry area
pN2O = 5*10^(-7); %mg/L %500ppb (0.5ppm) for N2O
PartialPressure =[pO2,pCO2,0,0,0,pNH3,pHONO,pN2O,pN2];
%tic
%[sitesC2ini, sitesCaini, sitesCgini, sitesCini] = Main_single_noBio_Nitrate(PartialPressure,amountCations,nitrateF, avgT, pot1);
%toc
tic
[sitesC2ini2, sitesCaini2, sitesCgini2, sitesCini2,sitesH2ONOini] = Main_single_noBio_nitroacidium(PartialPressure,amountCations,nitrateF, avgT, pot1);
toc

