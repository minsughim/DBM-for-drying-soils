
%Begin from 6AM (Sun rises t = 0)
%Oxygen KH
KHO2 = 1.3*0.001;
enthalpyO2 = 1500;
%Carbondioxide KH
KHCO2 = 3.5*0.01;
enthalpyCO2 = 2400;
%Ammonia KH
KHNH3 = 61;
enthalpyNH3 = 4200;
%HONO KH
KHHONO = 50;
enthalpyHONO = 4900;
%N2O KH
KHN2O = 2.5*0.01;
enthalpyN2O = 2600;

%HONO dissociation
KaHONO = 5.6*10^(-4);
%enthalpyHONOKa = 154.1; %calcuated myself
enthalpyHONOKa = 1423.5;
%enthalpyHONOKa = 715.39;
enthalpyHONOoverR = enthalpyHONOKa*1000/8.314;


stdTemp = 298.15;
absoluteT = 273.15;
conversionC = 12.2; % change the dimension to dimensionless henry's constant between aqueous phase and gas phase

pa = 1.01325e5;
%pa = 1e5;
partialPressureCO2 = 0.000383;
%partialPressureCO2 = 0.000500;
partialPressureO2 = 0.2095;


tempT = 100;
temperatureList = linspace(0,50,tempT);

for i = 1:tempT
        temperatureTemp = temperatureList(i);
        density_air = air_density(temperatureTemp, 0,pa);
        HenryConstListO2(i) = (temperatureTemp+absoluteT)*KHO2*exp(enthalpyO2*((1/(temperatureTemp+absoluteT))-(1/stdTemp)))/conversionC;
        HenryConstListCO2(i) = (temperatureTemp+absoluteT)*KHCO2*exp(enthalpyCO2*((1/(temperatureTemp+absoluteT))-(1/stdTemp)))/conversionC;
      %  HenryConstListNH3(i) = (temperatureTemp+absoluteT)*KHNH3*exp(enthalpyNH3*((1/(T))-(1/stdTemp)))/conversionC;
      %  HenryConstListHONO(i) = (temperatureTemp+absoluteT)*KHHONO*exp(enthalpyHONO*((1/(T))-(1/stdTemp)))/conversionC;
      %  HenryConstListN2O(i) = (temperatureTemp+absoluteT)*KHN2O*exp(enthalpyN2O*((1/(T))-(1/stdTemp)))/conversionC;

       O2solubility(i,1) = partialPressureO2*density_air*1000*HenryConstListO2(i);
       O2solubility(i,2) = partialPressureO2*1183.9*HenryConstListO2(i);      
       CO2solubility(i,1) = partialPressureCO2*density_air*1000*HenryConstListCO2(i);
       CO2solubility(i,2) = partialPressureCO2*1183.9*HenryConstListCO2(i);
       
       O2gas(i,1) = partialPressureO2*density_air*1000;
       O2gas(i,2) = partialPressureO2*1183.9;       
       CO2gas(i,1) = partialPressureCO2*density_air*1000;
       CO2gas(i,2) = partialPressureCO2*1183.9;
       
     %  CO2solubility(i,3) = partialPressureCO2*1215*HenryConstListCO2(i);
       
        
end