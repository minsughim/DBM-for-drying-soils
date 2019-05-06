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
Ynh4 = a(5); 
O2 = -1*a(2)*Molwt(2)/Ynh4; %O2
NH4 = -1*a(3)*Molwt(11)/Ynh4; %NH4+
CO2 = a(6)*Molwt(7)/Ynh4; %CO2
H2O = a(7)*Molwt(8)/Ynh4; %H2O
CH2O = -1*a(1)*Molwt(1)/Ynh4;
HCO3 = -1*a(4)*Molwt(10)/Ynh4;
NO3 = 0;
N2 = 0;
H = 0;
NO2 = 0;
N2O = 0;
%Yrespi = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);
YrespiNH4 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);
YrespiNH4_mole = [-1*a(2) a(6) -1*a(4) -1*a(1) NO3 -1*a(3) NO2 N2O N2 a(7) H]/Ynh4;


%% Phototrophs dark respiration  
%  
% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> O2 -> H2O : 
% + 0.25 O2 + 1 H+ + 1 e- --> + 0.5 H2O [ ?G = -78.72 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NO3-
% 0.1739 CO2 + 0.033 NO3- + 1.033 H+ + 1 e- --> 0.1739 C'1'H'1.6'O'0.4'N'0.19 + 0.3774 H2O [ ?G = 13.1803 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 O2 --> + 0.25 CO2 + 0.25 H2O  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.033 NO3- + 0.033 H+ --> 0.1739 C'1'H'1.6'O'0.4'N'0.19 + 0.0761 CO2 + 0.1274 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.26 ] [ fs = 0.74 ] [ e = 0.5 ]
% + 0.0417 C6H12O6 + 0.0651 O2 + 0.0244 NO3- + 0.0244 H+ --> 0.1286 C'1'H'1.6'O'0.4'N'0.19 + 0.1214 CO2 + 0.1593 H2O  
%  
% Yield prediction :
% Yg/m = 70.027 [ grams_cells/mol_donor ]
% Yc/m = 3.087 [ mol_C_cells/mol_donor ]
% Yc/c = 0.515 [ mol_C_cell/mol_C_donor ]


a = [0.0417*6 0.0651 0.0244 0.0244 0.1286 0.1214  0.1593];
Yno3 = a(5); 
O2 = -1*a(2)*Molwt(2)/Yno3; %O2
NH4 = 0; %NH4+
CO2 = a(6)*Molwt(7)/Yno3; %CO2
H2O = a(7)*Molwt(8)/Yno3; %H2O
CH2O = -1*a(1)*Molwt(1)/Yno3;
HCO3 = 0;
NO3 = -1*a(3)*Molwt(4)/Yno3;
N2 = 0;
H = -1*a(4)/Yno3;
NO2 = 0;
N2O = 0;
%Yrespi = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);
YrespiNO3 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);
YrespiNO3_mole = [-1*a(2) a(6) HCO3 -1*a(1) -1*a(3) NH4 NO2 N2O N2 a(7) -1*a(4)]/Yno3;

Y = Molwt(2)/Molwt(9);
%alphaCO = 4; %yield of C respected to the oxygen 
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
YH2Orequired = (YrespiNO3_mole(10) + alphaCO); % sugar required for respiration;
H2O = -1*YH2Orequired*Molwt(8)/Y;
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
YH2Orequired = (YrespiNH4_mole(10) + alphaCO); % sugar required for respiration;
H2O = -1*YH2Orequired*Molwt(8)/Y;
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
YH2Orequired = (YrespiNO3_mole(10) + alphaCO); % sugar required for respiration;
H2O = -1*YH2Orequired*Molwt(8)/Y;
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
YH2Orequired = (YrespiNH4_mole(10) + alphaCO); % sugar required for respiration;
H2O = -1*YH2Orequired*Molwt(8)/Y;
Yph4 = [O2 CO2 HCO3 CH2O NO3 NH4 NO2 N2O N2 H2O H]/Molwt(9);

Ymax(:,1) = Yph1(:);
Ymax(:,2) = Yph2(:);
Ymax(:,3) = Yph3(:);
Ymax(:,4) = Yph4(:);
