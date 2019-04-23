 function dy = HONO_forced_kinetics_Tdep_calcite(y,Z,pKs, theta)

% Y1 = NH4+
% Y2 = NH3
% Y3 = CO2
% Y4 = HCO3-
% Y5 = CO3 -2
% Y6 = NO2-
% Y7 = HONO
% Y8 = CaCO3 (aq)
% Y9 = CaCO3 (s)
% Y10 = Ca 2+
% Y11 = H +

%% Dynamics
dy(10,1) = 0;%zeros(10,1);    % a column vector
% Algebraic Equations : ionic interaction/tempearature dependency is not included in this case
Kaw = 10^(-14);
%Ka = 10^-9.4003;
%K1c = 10^-6.3819;
%K2c = 10^-10.3767;
%y(y<0) = 10^(-12);

y = real(y);

ions =[1 0 0 -1 -2 -1 0 0 0 2];
totIon = ions*y+Z; %Other cation and anoion included here for charge balancing
y11 = 0.5*(sqrt(totIon*totIon+4*Kaw)-totIon);

if y11 < Kaw %apporoximation of the concentraiton of OH-
    frac = sqrt(totIon*totIon+4*Kaw);
else
    frac = Kaw/y11;
end

I = 0.5*((ions.*ions)*y+Z+y11+frac); %ionic strength

absT = 273.15 + theta;
epsilon = 87.74-0.40008*theta + 9.398E-4*theta^2 -1.41E-6*theta^3;
B = 50.3/sqrt(epsilon*absT);
G = 1.825E+6*(epsilon*absT)^(-1.5);
% NH4+ HCO3- CO3^2- H+ Ca2+ NO2-
a = [2.5 4 4.5 9 6 3];
ionicValency = [1 -1 -2 1 2 -1];
%extended Debye-Heckel
pGamma = G.*ionicValency.^2.*(sqrt(I)./(1+ a.*B.*sqrt(I)));
%pGamma =zeros(6,1);

Ka = 10^(-pKs(1)+pGamma(4)-pGamma(1));
K1c = 10^(-pKs(2)+pGamma(4)+pGamma(2));
K2c = 10^(-pKs(3)+pGamma(4)+pGamma(3)-pGamma(2));
KcaComp = 10^(-pKs(4)+pGamma(3)+pGamma(5));
KcaPrecip = 10^(-pKs(5)+pGamma(3)+pGamma(5));
Kahono = 10^(-pKs(6)+pGamma(6)+pGamma(4));
%Kahono = 10^(-pKs(6)-pGamma(6)-pGamma(4));


%rate
k1 = 2221/(24*60*60); %/s
reac1 = k1*(y(3) - y(4)*y11/K1c); %CO2 + H2O = HCO3- + H+

%k2 = 7.19*10^8/(24*60*60); %/s
%pKa2 = 7.64; %very fast reaction -> needs to find tempearature dependency
%Ka2 = 10^(-pKs(3)+pGamma(4)+pGamma(3)-pGamma(2));
%reac2 = k2*(y(3)*frac-y(4)/Ka2); %CO2 + OH- = HCO3-
reac2 = 0; %ignore hydrolysis on OH-

%reacCO2 = reac1*(y(10)>K2c) +reac2*(y(10)>K2c);
k3 = 10^10/(24*60*60); %/s
reac3 = k3*(y(4)-y11*y(5)/K2c);%HCO3- = H+ +CO32-
k4 = 10^10/(24*60*60); %/s
reac4 = k4*(y(1)-y11*y(2)/Ka);%NH4 = NH3 + H+

% Nitrous acid and nitrite reaction
k5 = 10^10/(24*60*60); %/s % Horokawa says 10^(6)s^(-1) almost 100 times faster than the value used here
reac5 = k5*(y(7)-y11*y(6)/Kahono);%HNO2 = NO2- + H+

% complexation of calcite
kcom = 10^10/(24*60*60); %/s assuming fast reaction
reac6 = kcom*(y(8)-KcaComp*y(5)*y(10));

% Precipitayion of calcite
kpre = 10^10/(24*60*60); %/s assuming fast reaction
reac7 = kpre*(y(9)-y(5)*y(10)/KcaPrecip);


% Y5 = CO3 -2
% Y6 = NO2-
% Y7 = HONO
% Y8 = CaCO3 (aq)
% Y9 = CaCO3 (s)
% Y10 = Ca

dy(1) = -reac4;
dy(2) = reac4;
dy(3) = -reac1-reac2;
dy(4) = reac1+reac2-reac3;
dy(5) = reac3+reac6+reac7;
dy(6) = reac5;
dy(7) = -reac5;
dy(8) = -reac6;
dy(9) = -reac7;
dy(10) = reac6+reac7;

% if (y(3) > kH && dy(3)>0)
%    dy(3) = 0; 
%    dy(4) = -reac3;
% end


end