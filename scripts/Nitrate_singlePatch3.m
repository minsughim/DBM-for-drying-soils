clear all
maxT = 2;
poolobj = parpool('local',maxT);

pCO2frac = 1;
amountCations = 0.1;
avgT = 25;
HONOppb = 0.005;
NH3ppb = 20;
pN2 = 0.7809;
pO2 = 0.2095;
pCO2 = 0.000383*pCO2frac; %383 ppm
pNH3 = NH3ppb*10^(-9); %default: 5ppb: Gong et al 2011 %2ppb in McCally
pHONO = HONOppb*10^(-9); %1ppb for HONO : Su et al 2011 from dry area
pN2O = 5*10^(-7); %mg/L %500ppb (0.5ppm) for N2O
PartialPressure =[pO2,pCO2,0,0,0,pNH3,pHONO,pN2O,pN2];
potResolution = 5;
nitrateResolution = 5;

potList = logspace(0,3,potResolution);
nitrateList = linspace(0,10,nitrateResolution);
pHList = zeros(nitrateResolution,potResolution);
WFList = zeros(nitrateResolution,potResolution);
GCList = zeros(nitrateResolution,potResolution);


for i = 1:nitrateResolution
    tic
    nitrateF = nitrateList(i);
    i
    parfor j = 1:potResolution
        pot1 = potList(j);
        [sitesC2ini, sitesCaini, sitesCgini, sitesCini, sitesH2ONOini, waterFilm, LocalGasContents] = Main_single_noBio_nitroacidium_test_Patch(PartialPressure,amountCations,nitrateF,avgT,pot1);
        pHList(i,j) = sitesC2ini{5};
        WFList(i,j) = waterFilm;
        GCList(i,j) = LocalGasContents;
    end
    toc
end

serial_id = sprintf('Single_pHdependence_NH3%d_HONO%d_Cat%.1f',NH3ppb, HONOppb,amountCations);
save(strcat(serial_id,'.mat'))

delete(poolobj)
