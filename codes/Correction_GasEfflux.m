examineDay = 5;
potS = 3;
plottt = 5;
maximumThreadt = maxT;
pCO2frac = 1;
amountCations = 0.1;
avgT = 25;
alphaCO = 3;
lightOn = 0;
%indexS = 7;


drypot2List = [50, 100, 500, 10000];

for indexPot = 1:4
    
    drypot2 = drypot2List(indexPot);

    serial_id_new = sprintf('HT_BSC_photoY3_Day%dPot%.1fpCO2frac%dCat%.1favT%dAlphaCO%dfindex%d',examineDay,potS,pCO2frac,amountCations,avgT,alphaCO, indexS);
cd(serial_id_new);
         serial_idtemp = sprintf('Constant_dark_water_decay4_Drying%dkPa_varyingT', drypot2);
   
   cd(serial_idtemp); 
  load('Final.mat')
  cd ..
  cd ..
  
  EffluxCalculation_dynamics
  totVerticalArea = mean(min(totaldM))*Lp*n;
  GasEffluxResult = effluxList2*10^6/dt/plottt/totVerticalArea; %\mug/m^2s
  save(strcat('./',serial_id_new,'/',serial_idtemp,'/GasEfflux.mat'),'effluxList2', 'GasEffluxResult','totVerticalArea','timeSeriesDeltaM');
  
end

