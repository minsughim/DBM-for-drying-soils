numberOFN = 8; %O2, CO2, HCO3-, sugar, NO3-, NH4+, NO2- N2O: used for food or byproduct
numberOFP = 8; % four phototrphs(nitrogen fixing, subindexed), one aerobe, one anaerobe (denitrifying), one chemoautotrophic AOB, and NOB.
%Chemical Factor (food for microbes)
%molecular weight information for the stoichiometry
% Species                    Mol. Wt.
% -------                    --------
% 1: CH2O                          30.03
% 2: O2                            32.00
% 3: NH3                           17.03
% 4: NO3                           62.00
% 5: N2                            28.01
% 6: CH1.8O0.5N0.2                 24.63
% 7: CO2                           44.01
% 8: H2O                           18.02
% 9: CH2.5ON0.17                 32.91
%10: HCO3                          61.02
%11: NH4                           18.04
%12: CO3                           60.01
%13: NO2                           46
%14: N2O                           28.01+16
Molwt = [30.0263 31.9988 17.0306 62.0049 28.0135 24.6263 44.0098 18.0153 32.9114 61.0171 18.0385 60.0092 62 28.01+16];
%% nutrient description

%N1 = O2
%Carbon related
%N2 = CO2
%N3 = HCO3-
%N4 = carbohydrates (CH2O) : sugar like glucose
%Nitrogen related
%N5 = NO3-
%N6 = NH3
%N7 = N2

iniConcentrationGas = zeros(numberOFN,1);
iniConcentrationGas(1) =  sitesCgini{1}; %mg/L
iniConcentrationGas(2) =  sitesCgini{2}; %mg/L
iniConcentrationGas(6) =  sitesCgini{6}; %mg/L %atmospheric level of ammonia
iniConcentrationGas(7) =  sitesCgini{7}; %mg/L % for HONO
iniConcentrationGas(8) = sitesCgini{8}; %mg/L %100ppb for N2O

%% All dissolved
iniConcentration(1) = sitesCini{1}; % N4 : O2 mg/L : Henry's law at 1 atm
iniConcentration(2) = sitesCini{2}; % N1 : CO2 mg/L : Henry's law at 1 atm
iniConcentration(3) = sitesCini{3}; % at equilbrium with atmospheric level
iniConcentration(4) = sitesCini{4}; % N5 : SS 1M (readily degradable organic substrate) : CH1.5O0.5N0.1 %Initmole only controls the sugar concentraion. (Catillo-Monroy 2010): 34.61mg/kg soil: with porosity 0.3: 34.61*2.6*0.3~
iniConcentration(5) = sitesCini{5}; % N2 : NO3- : (Abed 2013, for hot desert with lots of nitrate)
iniConcentration(6) = sitesCini{6}; %start with charge neutral condition without buffer
iniConcentration(7) = sitesCini{7}; %NO2-
iniConcentration(8) = sitesCini{8}; %N2O

%iniConcentration(6) = 1.755*initMoles; % N3 : NH3 +NH4 :  total Ammonium: (Catillo-Monroy 2010): 2.25mg/kg soil: with porosity 0.3: 2.25*2.6*0.3~ 1.7
EPSCcrit= 2000; % EPS gellation at 0.2%
%EPSMcrit = EPSCcrit*10^6/(1-EPSCcrit);
% %%%%%%%%

%% Diffusion coefficient
%[O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O]
%Fix the diffusion coefficient for EPS later 
DiffusionList = [1.73 1.65 1.02 0.5184 1.47 1.41 1.64 1.58]*10^(-4)/(24*3600); %m^2/s %nitrate diffusion coefficient from Hirokawa
%[EPS CO3 NH3 H+ pH HONO Ca2+]
DiffusionAuxil = [0.1 0.79 1.41 8.04 0 1.64 0.864]*10^(-4)/(24*3600); %m^2/s  

%% Parameters for the microbial life
V0r = 1/14;
chi0r = 1/14;
V0 = V0r*0.000014;
chi0 = chi0r*50*1e-9;
Vu      =   10*0.4*1e-18;%median cell volume at mu=0, unit:[fl fl].*1e-9 = mm^3
R       =   1e-6; % radius of the bacteria : assume that its shape is cylindrical
rho     =   290*(10^3);%cell density (dry mass), unit:[fg/fl fg/fl].*1e3 = ng/mm^-3
Vdmin   =   Vu*(2/1.433);%minimum cell volume at division, unit:[fl fl].*1e-9 = mm^-3, Descriptive model
Vmin    =   Vdmin/5;%minimum cell volume, assumed to be 1/5 of Vdmin
germT   =   1/8; %days
mumax = 3*[10 10 10 10 5.5 1.6 1.7 1]/(24*60*60); %Maximum specific growth rate, unit: (s-1) % AOB from bollmann et al
mumax(end+1) = mumax(1)/10;
mrate   =  3*[0.09 0.09 0.09 0.09 0.4 0.4 0.15 0.15]/(24*60*60);%maintenance rate, unit :per second

pGmax = 1000;
pGmin = 0.01;

KsP = 1.4;
KsI = 295;

NfixationR = 0.1; % Nitrgoen fixation
repsiRatio = 0.1; % respiration rate of photoautotrophs is proportional to  
germT = 1; % 1 hour: response delay rather than the germination time


%% Stoichiometry
multiplytingF = 1; % In case to consider that Cyanobacteria has a higher carbon contents per cell compare to other bacterial cells 
Molwt(9) = multiplytingF*22.684437; %C'1'H'1.6'O'0.4'N'0.19
cph = multiplytingF*1;
hph = multiplytingF*1.6;
oph = multiplytingF*0.4;
nph = multiplytingF*0.19;


%% Nitrogen fixation by cyanobacteria
% %  
% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> N2 -> NH4+ : 
% + 0.1667 N2 + 1.3333 H+ + 1 e- --> + 0.3333 NH4+ [ ?G = 26.7 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.1915 CO2 + 0.0449 HCO3- + 0.0449 NH4+ + 1 H+ + 1 e- --> 0.2364 C'1'H'1.6'O'0.4'N'0.19 + 0.4232 H2O [ ?G = 17.8253 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 H2O + 0.1667 N2 + 0.3333 H+ --> + 0.25 CO2 + 0.3333 NH4+  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0449 NH4+ + 0.0449 HCO3- --> 0.2364 C'1'H'1.6'O'0.4'N'0.19 + 0.0585 CO2 + 0.1732 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.81 ] [ fs = 0.19 ] [ e = 0.5 ]
% + 0.0417 C6H12O6 + 0.0087 HCO3- + 0.1679 H2O + 0.1343 N2 + 0.2687 H+ --> 0.0459 C'1'H'1.6'O'0.4'N'0.19 + 0.2128 CO2 + 0.2599 NH4+  
%  
% Yield prediction :
% Yg/m = 24.973 [ grams_cells/mol_donor ]
% Yc/m = 1.101 [ mol_C_cells/mol_donor ]
% Yc/c = 0.183 [ mol_C_cell/mol_C_donor ]

a = [0.0417*6 0.0087 0.1679 0.1343 0.2687 0.0459 00.2128 0.2599];  
Y = a(6);
O2 = 0; %O2
NH4 = a(8)*Molwt(11)/Y; %NH4
CO2 = a(7)*Molwt(7)/Y; %CO2
H2O = -1*a(3)*Molwt(8)/Y; %H2O
CH2O = -1*a(1)*Molwt(1)/Y;
HCO3 = -1*a(2)*Molwt(10)/Y;
NO3 =  0;
N2 = -a(4)*Molwt(5)/Y;
H = -a(5)/Y;
NO2 = 0;
N2O =0;
YN2fix = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);


Y = Molwt(2)/Molwt(9);
alphaCO = 5; %yield of C respected to the oxygen for phototrophs
BiomassY = [cph+alphaCO; hph+2*alphaCO; oph+alphaCO+2; nph];

%a1CO2 +a2NO3 +a3H2O +a4H+ =YCH2.5ON0.17 + alphaCH2O+O2
temp = [1 0 0 0; 0 0 2 1; 2 3 1 0; 0 1 0 0];% CO2 NO3-
a = temp\BiomassY;
O2 = Molwt(2)/Y; %O2
NH4 = 0;  %NH4+
CO2 = -1*a(1)*Molwt(7)/Y; %CO2
H2O = -1*a(3)*Molwt(8)/Y; %H2O
CH2O = alphaCO*Molwt(1)/Y;
HCO3 = 0;
NO3 =  -1*a(2)*Molwt(4)/Y;
N2 = 0;
H = -1*a(4)/Y;
NO2 = 0;
N2O = 0;
Yph1 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);

%a1CO2 +a2NH4+ +a3H2O +a4H+ =YCH2.5ON0.17+ alphaCH2O+O2
temp = [1 0 0 0; 0 4 2 1; 2 0 1 0; 0 1 0 0]; %CO2 NH4+
a = temp\BiomassY;
O2 = Molwt(2)/Y; %O2
NH4 = -1*a(2)*Molwt(11)/Y;  %NH4+
CO2 = -1*a(1)*Molwt(7)/Y; %CO2
H2O = -1*a(3)*Molwt(8)/Y; %H2O
CH2O = alphaCO*Molwt(1)/Y;
HCO3 = 0;
NO3 =  0;
N2 = 0;
H = -1*a(4)/Y;
NO2 = 0;
N2O = 0;
Yph2 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);

%a1HCO3 +a2NO3 +a3H2O +a4H+ =YCH2.5ON0.17 + alphaCH2O+O2
temp = [1 0 0 0; 1 0 2 1; 3 3 1 0; 0 1 0 0]; %HCO3- NO3-
a = temp\BiomassY;
O2 = Molwt(2)/Y; %O2
NH4 = 0;  %NH4+
CO2 = 0; %CO2
H2O = -1*a(3)*Molwt(8)/Y; %H2O
CH2O = alphaCO*Molwt(1)/Y;
HCO3 = -1*a(1)*Molwt(10)/Y;
NO3 =  -1*a(2)*Molwt(4)/Y;
N2 = 0;
H = -1*a(4)/Y;
NO2 = 0;
N2O = 0;
Yph3 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);

%a1HCO3- +a2NH4+ +a3H2O +a4H+ =YCH2.5ON0.17+ alphaCH2O+O2
temp = [1 0 0 0; 1 4 2 1; 3 0 1 0; 0 1 0 0]; %HCO3- NH4+
a = temp\BiomassY;
O2 = Molwt(2)/Y; %O2
NH4 = -1*a(2)*Molwt(11)/Y;  %NH4+
CO2 = 0; %CO2
H2O = -1*a(3)*Molwt(8)/Y; %H2O
CH2O = alphaCO*Molwt(1)/Y;
HCO3 = -1*a(1)*Molwt(10)/Y;
NO3 =  0;
N2 = 0;
H = -1*a(4)/Y;
NO2 = 0;
N2O = 0;
Yph4 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);

%% Phototrophs dark respiration  

% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> O2 -> H2O : 
% + 0.25 O2 + 1 H+ + 1 e- --> + 0.5 H2O [ ?G = -78.72 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.1915 CO2 + 0.0449 HCO3- + 0.0449 NH4+ + 1 H+ + 1 e- --> 0.2364 C'1'H'1.6'O'0.4'N'0.19 + 0.4232 H2O [ ?G = 17.8253 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 O2 --> + 0.25 CO2 + 0.25 H2O  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0449 NH4+ + 0.0449 HCO3- --> 0.2364 C'1'H'1.6'O'0.4'N'0.19 + 0.0585 CO2 + 0.1732 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.34 ] [ fs = 0.66 ] [ e = 0.5 ]
% + 0.0417 C6H12O6 + 0.0841 O2 + 0.0298 NH4+ + 0.0298 HCO3- --> 0.1569 C'1'H'1.6'O'0.4'N'0.19 + 0.1229 CO2 + 0.199 H2O  
%  
% Yield prediction :
% Yg/m = 85.415 [ grams_cells/mol_donor ]
% Yc/m = 3.765 [ mol_C_cells/mol_donor ]
% Yc/c = 0.628 [ mol_C_cell/mol_C_donor ]

a = [0.0417*6 0.0841 0.0298 0.0298 0.1569 0.1229 0.199];
Y = a(5); 
O2 = -1*a(2)*Molwt(2)/Y; %O2
NH4 = -1*a(3)*Molwt(11)/Y; %NH4+
CO2 = a(6)*Molwt(7)/Y; %CO2
H2O = a(7)*Molwt(8)/Y; %H2O
CH2O = -1*a(1)*Molwt(1)/Y;
HCO3 = -1*a(4)*Molwt(10)/Y;
NO3 = 0;
N2 = 0;
H = 0;
NO2 = 0;
N2O = 0;
Yrespi = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);

%% Heterotrohps(aerobes)  

% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> O2 -> H2O : 
% + 0.25 O2 + 1 H+ + 1 e- --> + 0.5 H2O [ ?G = -78.72 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.1905 CO2 + 0.0476 HCO3- + 0.0476 NH4+ + 1 H+ + 1 e- --> 0.2381 C'1'H'1.8'O'0.5'N'0.2 + 0.4048 H2O [ ?G = 19.4892 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 O2 --> + 0.25 CO2 + 0.25 H2O  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0476 NH4+ + 0.0476 HCO3- --> 0.2381 C'1'H'1.8'O'0.5'N'0.2 + 0.0595 CO2 + 0.1548 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.36 ] [ fs = 0.64 ] [ e = 0.5 ]
% + 0.0417 C6H12O6 + 0.09 O2 + 0.0305 NH4+ + 0.0305 HCO3- --> 0.1524 C'1'H'1.8'O'0.5'N'0.2 + 0.1281 CO2 + 0.189 H2O  
%  
% Yield prediction :
% Yg/m = 90.074 [ grams_cells/mol_donor ]
% Yc/m = 3.658 [ mol_C_cells/mol_donor ]
% Yc/c = 0.61 [ mol_C_cell/mol_C_donor ]

a = [0.0417*6 0.09 0.0305  0.0305 0.1524 0.1281 0.189];
Y = a(5); 
O2 = -1*a(2)*Molwt(2)/Y; %O2
NH4 = -1*a(3)*Molwt(11)/Y; %NH4+
CO2 = a(6)*Molwt(7)/Y; %CO2
H2O = a(7)*Molwt(8)/Y; %H2O
CH2O = -1*a(1)*Molwt(1)/Y;
HCO3 = -1*a(4)*Molwt(10)/Y;
NO3 = 0;
N2 = 0;
H = 0;
NO2 = 0;
N2O = 0;
Yaerobe = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(6);

%% Heterotrohps(anaerobes: denitrification)  

% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> NO3- -> N2O : 
% + 0.25 NO3- + 1.25 H+ + 1 e- --> + 0.125 N2O + 0.625 H2O [ ?G = -57.54 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NO3-
% 0.1724 CO2 + 0.0345 NO3- + 1.0345 H+ + 1 e- --> 0.1724 C'1'H'1.8'O'0.5'N'0.2 + 0.3621 H2O [ ?G = 14.1851 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 NO3- + 0.25 H+ --> + 0.25 CO2 + 0.125 N2O + 0.375 H2O  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0345 NO3- + 0.0345 H+ --> 0.1724 C'1'H'1.8'O'0.5'N'0.2 + 0.0776 CO2 + 0.1121 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.32 ] [ fs = 0.68 ] [ e = 0.5 ]
% + 0.0417 C6H12O6 + 0.1032 NO3- + 0.1032 H+ --> 0.1174 C'1'H'1.8'O'0.5'N'0.2 + 0.1326 CO2 + 0.0399 N2O + 0.1959 H2O  
%  
% Yield prediction :
% Yg/m = 69.408 [ grams_cells/mol_donor ]
% Yc/m = 2.819 [ mol_C_cells/mol_donor ]
% Yc/c = 0.47 [ mol_C_cell/mol_C_donor ]
a = [0.0417*6 0.1032 0.1032 0.1174 0.1326 0.0399 0.1959];
Y = a(4); 
CH2O = -a(1)*Molwt(1)/Y;
NO3  = -1*a(2)*Molwt(4)/Y;
CO2 = a(5)*Molwt(7)/Y;
N2 = 0;
H2O = a(7)*Molwt(8)/Y;
O2 = 0;
HCO3 = 0;
NH4 = 0;
H = -1*a(3)/Y;
NO2 = 0;
N2O = a(6)*Molwt(14)/Y;
Ydenitri = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(6);


%% Chemoautotrophs(nitrifier1): convert ammonium to nitrite (AOB)
% (rd) Electron donor --> NH4+ -> NO2- : 
% + 0.1667 NO2- + 1.3333 H+ + 1 e- --> + 0.1667 NH4+ + 0.3333 H2O [ ?G = -32.93 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> O2 -> H2O : 
% + 0.25 O2 + 1 H+ + 1 e- --> + 0.5 H2O [ ?G = -78.72 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.1905 CO2 + 0.0476 HCO3- + 0.0476 NH4+ + 1 H+ + 1 e- --> 0.2381 C'1'H'1.8'O'0.5'N'0.2 + 0.4048 H2O [ ?G = 19.4892 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.25 O2 + 0.1667 NH4+ --> + 0.3333 H+ + 0.1667 NO2- + 0.1667 H2O  
%  
% Synthesis reaction : 
% + 0.1905 CO2 + 0.2143 NH4+ + 0.0476 HCO3- --> 0.2381 C'1'H'1.8'O'0.5'N'0.2 + 0.3333 H+ + 0.1667 NO2- + 0.0714 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.88 ] [ fs = 0.12 ] [ e = 0.5 ]
% + 0.2198 O2 + 0.023 CO2 + 0.1724 NH4+ + 0.0058 HCO3- --> 0.0288 C'1'H'1.8'O'0.5'N'0.2 + 0.3333 H+ + 0 e- + 0.1667 NO2- + 0.1552 H2O  
%  
% Yield prediction :
% Yg/m = 4.108 [ grams_cells/mol_donor ]
% Yc/m = 0.167 [ mol_C_cells/mol_donor ]

a = [0.2198 0.023 0.1724 0.0058 0.0288 0.3333 0.1667 0.1552];
Y = a(5); 
O2 = -1*a(1)*Molwt(2)/Y;
CO2 = -1*a(2)*Molwt(7)/Y;
NH4 = -1*a(3)*Molwt(11)/Y;
NO3 = 0;
H2O = a(8)*Molwt(8)/Y;
HCO3 = -1*a(4)*Molwt(10)/Y;
CH2O = 0;
N2 = 0;
H = a(6)/Y;
NO2 = a(7)*Molwt(13)/Y;
N2O = 0;
Ynitri1 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(6);

%% Chemoautotrophs(nitrifier): nitrite oxidizing bacteria (NOB)

% (rd) Electron donor --> NO2- -> NO3- : 
% + 0.5 NO3- + 1 H+ + 1 e- --> + 0.5 NO2- + 0.5 H2O [ ?G = -41.65 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> O2 -> H2O : 
% + 0.25 O2 + 1 H+ + 1 e- --> + 0.5 H2O [ ?G = -78.72 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NO2-
% 0.1852 CO2 + 0.037 NO2- + 1.037 H+ + 1 e- --> 0.1852 C'1'H'1.8'O'0.5'N'0.2 + 0.3519 H2O [ ?G = 15.1993 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.25 O2 + 0.5 NO2- --> + 0.5 NO3-  
%  
% Synthesis reaction : 
% + 0.1852 CO2 + 0.1481 H2O + 0.537 NO2- + 0.037 H+ --> 0.1852 C'1'H'1.8'O'0.5'N'0.2 + 0.5 NO3-  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.9 ] [ fs = 0.1 ] [ e = 0.5 ]
% + 0.2261 O2 + 0.0177 CO2 + 0.0142 H2O + 0.5035 NO2- + 0.0035 H+ --> 0.0177 C'1'H'1.8'O'0.5'N'0.2 + 0.5 NO3-  
%  
% Yield prediction :
% Yg/m = 0.865 [ grams_cells/mol_donor ]
% Yc/m = 0.035 [ mol_C_cells/mol_donor ]

a = [0.2261 0.0177 0.0142 0.5035 0.0035 0.0177 0.5];
Y = a(6); 
O2 = -a(1)*Molwt(2)/Y;
CO2 = -1*a(2)*Molwt(7)/Y;
NH4 = 0;
NO3 = a(7)*Molwt(4)/Y;
H2O = -1*a(3)*Molwt(8)/Y;
HCO3 = 0;
CH2O = 0;
N2 = 0;
H = -1*a(5)/Y;
NO2 = -1*a(4)*Molwt(13)/Y;
N2O = 0;
Ynitri2 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(6);


%%
Y = 0.5; %yield of 50% going back to the environment
YDecay = zeros(2,11);
%%%%%
%decay of autotrophs
x = 1.6;y = 0.4;z=0.19;
A =[0 0 1 0 1;-2 4 1 1 x; 1 0 -3 0 -y;0 1 0 0 z; 0 1 -1 1 0];
invA = inv(A);
temp = [1-Y; x-2*Y; Y-y; z; 0];
a = invA*temp;
%CHxOyNz +a1H2O = YCH2O + a2NH4+ + a3HCO3- + a4H+ +a5CHxOyNz
NH4 = a(2)*Molwt(11)/a(5);
H2O = -1*a(1)*Molwt(8)/a(5);
HCO3 = a(3)*Molwt(10)/a(5);
CH2O = Y/a(5);
H = a(4)/a(5);
CHxOyNz = Molwt(9);
YDecay(1,:) = [0 0 HCO3 CH2O 0 NH4 0 0 0 H2O H]/Molwt(9);

%decay of heterotrophs
x = 1.8;y = 0.5;z=0.2;
A =[0 0 1 0 1;-2 4 1 1 x; 1 0 -3 0 -y;0 1 0 0 z; 0 1 -1 1 0];
invA = inv(A);
temp = [1-Y; x-2*Y; Y-y; z; 0];
a = invA*temp;
%CHxOyNz +a1H2O = YCH2O + a2NH4+ + a3HCO3- + a4H+ +a5CHxOyNz
NH4 = a(2)*Molwt(11)/a(5);
H2O = -1*a(1)*Molwt(8)/a(5);
HCO3 = a(3)*Molwt(10)/a(5);
CH2O = Y/a(5);
H = a(4)/a(5);
CHxOyNz = Molwt(6);
YDecay(2,:) = [0 0 HCO3 CH2O 0 NH4 0 0 0 H2O H]/Molwt(6);

Ymax = zeros(numberOFP, numberOFN+3);
Ymax(1,:) = Yph1;
Ymax(2,:) = Yph2;
Ymax(3,:) = Yph3;
Ymax(4,:) = Yph4;
Ymax(5,:) = Yaerobe;
Ymax(6,:) = Ydenitri;
Ymax(7,:) = Ynitri1;
Ymax(8,:) = Ynitri2;

YmaxCNratio(1) = Ymax(1,2)/Ymax(1,5);
YmaxCNratio(2) = Ymax(2,2)/Ymax(2,6);
YmaxCNratio(3) = Ymax(3,3)/Ymax(3,5);
YmaxCNratio(4) = Ymax(4,3)/Ymax(4,6);

Ymax = Ymax';
YmaxCandN = Ymax;

EPSfraction = 0.2;
YrespiTot = [Yrespi; Yrespi; Yrespi; Yrespi; 0*Yrespi; 0*Yrespi; 0*Yrespi; 0*Yrespi]';


%% Growth rate with Ks values 

PhotoRespKO = 9.599; % amount of oxygen outside oc the cell
PhotoRespSugar = 3.411285;
%Oxygen half-saturation is fixed based on Henze et al 2000 -> activated
%sludge model. (Wolf 2007 has a typo.. it seems.
oxygenK = 5*0.2496;

%[O2 CO2 HCO3 CH2O NO3 NH3 NO2] %HCO3 for KS is not given in wolf et al.
%assumed that it is the same as CO2 Km
KsPh1 = [PhotoRespKO 0.88 0 PhotoRespSugar 0.0744 -0.02 0];
KsPh2 = [PhotoRespKO 0.88 0 PhotoRespSugar 0 0.02 0];
KsPh3 = [PhotoRespKO -0.0088 0.88 PhotoRespSugar 0.0744 -0.02 0];
KsPh4 = [PhotoRespKO -0.0088 0.88 PhotoRespSugar 0 0.02 0];
HeteroKs = 10^(-3);


Ks = [KsPh1;KsPh2;KsPh3;KsPh4;oxygenK 0 0.88 HeteroKs 0 1.7*10^(-6) 0;-1*oxygenK 0 0 HeteroKs 0.0022 0 0; 0.5 0.88 0 0 0 0.9 0; 0.5 0.88 0 0 0 0 2.96];%Nitrobactor's Ks value for NO2- has been used, Nitrosonma's minimum ammonium Km

