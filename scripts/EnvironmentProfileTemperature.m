function [intensityList,MumaxT, HenryConstListO2,HenryConstListCO2,HenryConstListNH3, HenryConstListHONO, HenryConstListN2O, pKList,Density_air] = EnvironmentProfileTemperature(depthList,timeLine,sitesTevolove)

%Calculate Henry's constant (Added 11th March 2014)

global attenuationRate MaximumIncidence m n

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
enthalpyHONOKa = 1423.5;
enthalpyHONOoverR = enthalpyHONOKa*1000/8.314;


stdTemp = 298.15;
absoluteT = 273.15;
conversionC = 12.2; % change the dimension to dimensionless henry's constant between aqueous phase and gas phase

tempT = length(timeLine);
tempL =length(depthList);
intensityList = zeros(tempT, tempL);
MumaxT = zeros(m,n,tempT);
HenryConstListO2 = zeros(m,n,tempT);
HenryConstListCO2 = zeros(m,n,tempT);
HenryConstListNH3 = zeros(m,n,tempT);
HenryConstListHONO = zeros(m,n,tempT);
HenryConstListN2O = zeros(m,n,tempT);
Density_air = zeros(m,n,tempT);
pKList = cell(m,n,tempT); %Ka K1c, K2c + calcite precipitation, association + HONO and nitrite
%Calcuate temperatue depedent mumax
%b = 0.0410;
%Tmin = 4;
%c = 0.161;
%Tmax = 43.7;

%mu25 = 1.42;
Tl = 297.7;
Th = 314.7;
Hh = 687900;
Hl = -141100;
Ha = -5430;

%Mum = (mu25/stdTemp)*T*exp((Ha/8.31)*(1/stdTemp -1/T))/(1+exp((Hl/8.31)*(1/Tl -(1/T)))+exp((Hh/8.31)*(1/Th -(1/T))))
pa = 1.01325e5; % 1 atm

for timeT = 1:tempT
    
    tTemp = timeLine(timeT);
    
    for i = 1:m
        
        zTemp = depthList(i);
        
        for j = 1:n
            
            temperatureTemp = sitesTevolove(i,j,timeT);
            T = temperatureTemp+absoluteT;
            pKList{i,j,timeT}(1) = 2835.76/T - 0.6322 + 0.001225*T;
            pKList{i,j,timeT}(2) = 3404.71/T - 14.8435 + 0.032786*T;
            pKList{i,j,timeT}(3) = 2902.39/T - 6.498 + 0.02379*T;
            pKList{i,j,timeT}(4) = 171.9065 + 0.077993*T - 2839.319/T - 71.595*log10(T);
            pKList{i,j,timeT}(5) = 1228.732 + 0.299444*T - 35512.75/T - 485.818*log10(T);
            pKList{i,j,timeT}(6) = -log10(KaHONO*exp(enthalpyHONOoverR*((1/T)-(1/stdTemp))));
            
            MumaxT(i,j,timeT) = (1/stdTemp)*T*exp((Ha/8.31)*(1/stdTemp -1/T))/(1+exp((Hl/8.31)*(1/Tl -(1/T)))+exp((Hh/8.31)*(1/Th -(1/T))));
            HenryConstListO2(i,j,timeT) = (temperatureTemp+absoluteT)*KHO2*exp(enthalpyO2*((1/(T))-(1/stdTemp)))/conversionC;
            HenryConstListCO2(i,j,timeT) = (temperatureTemp+absoluteT)*KHCO2*exp(enthalpyCO2*((1/(T))-(1/stdTemp)))/conversionC;
            HenryConstListNH3(i,j,timeT) = (temperatureTemp+absoluteT)*KHNH3*exp(enthalpyNH3*((1/(T))-(1/stdTemp)))/conversionC;
            HenryConstListHONO(i,j,timeT) = (temperatureTemp+absoluteT)*KHHONO*exp(enthalpyHONO*((1/(T))-(1/stdTemp)))/conversionC;
            HenryConstListN2O(i,j,timeT) = (temperatureTemp+absoluteT)*KHN2O*exp(enthalpyN2O*((1/(T))-(1/stdTemp)))/conversionC;
            
            Density_air(i,j,timeT) = air_density(temperatureTemp, 0,pa)*1000; %mg/L
        end
        
        if (mod(tTemp,24)<12)&&(mod(tTemp,24)~=0)
            intensityList(timeT,i) = (0.5*(1-cos(2*pi*tTemp/12)))*MaximumIncidence.*exp(-zTemp/(attenuationRate));
        else
            intensityList(timeT,i) = 0;
        end
        
        
    end
end


end