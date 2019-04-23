% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> N2 -> NH4+ : 
% + 0.1667 N2 + 1.3333 H+ + 1 e- --> + 0.3333 NH4+ [ ?G = 26.7 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.208 CO2 + 0.0426 HCO3- + 0.0426 NH4+ + 1 H+ + 1 e- --> 0.2506 C'1'H'2.5'O'1'N'0.17 + 0.2932 H2O [ ?G = 27.4176 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 H2O + 0.1667 N2 + 0.3333 H+ --> + 0.25 CO2 + 0.3333 NH4+  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0426 NH4+ + 0.0426 HCO3- --> 0.2506 C'1'H'2.5'O'1'N'0.17 + 0.042 CO2 + 0.0432 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.92 ] [ fs = 0.08 ] [ e = 0.4 ]
% + 0.0417 C6H12O6 + 0.0036 HCO3- + 0.2255 H2O + 0.1528 N2 + 0.3055 H+ --> 0.0209 C'1'H'2.5'O'1'N'0.17 + 0.2326 CO2 + 0.302 NH4+  
%  
% Yield prediction :
% Yg/m = 16.519 [ grams_cells/mol_donor ]
% Yc/m = 0.502 [ mol_C_cells/mol_donor ]
% Yc/c = 0.084 [ mol_C_cell/mol_C_donor ]

a = [0.0417*6 0.0036 0.2255 0.1528 0.3055 0.0209 0.2326 0.302];  
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
% 
% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> N2 -> NH4+ : 
% + 0.1667 N2 + 1.3333 H+ + 1 e- --> + 0.3333 NH4+ [ ?G = 26.7 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.208 CO2 + 0.0426 HCO3- + 0.0426 NH4+ + 1 H+ + 1 e- --> 0.2506 C'1'H'2.5'O'1'N'0.17 + 0.2932 H2O [ ?G = 27.4176 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 H2O + 0.1667 N2 + 0.3333 H+ --> + 0.25 CO2 + 0.3333 NH4+  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0426 NH4+ + 0.0426 HCO3- --> 0.2506 C'1'H'2.5'O'1'N'0.17 + 0.042 CO2 + 0.0432 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.82 ] [ fs = 0.18 ] [ e = 0.6 ]
% + 0.0417 C6H12O6 + 0.0078 HCO3- + 0.1965 H2O + 0.1363 N2 + 0.2726 H+ --> 0.0457 C'1'H'2.5'O'1'N'0.17 + 0.2121 CO2 + 0.2648 NH4+  
%  
% Yield prediction :
% Yg/m = 36.089 [ grams_cells/mol_donor ]
% Yc/m = 1.097 [ mol_C_cells/mol_donor ]
% Yc/c = 0.183 [ mol_C_cell/mol_C_donor ]


a = [0.0417*6 0.0078  0.1965 0.1363 0.2726 0.0457 0.2121 0.2648];  
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


 
% (rd) Electron donor --> Glucose : 
% + 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ?G = 41.35 KJ/e-eq ] 
%  
% (ra) Electron acceptor --> N2 -> NH4+ : 
% + 0.1667 N2 + 1.3333 H+ + 1 e- --> + 0.3333 NH4+ [ ?G = 26.7 KJ/e-eq ] 
%  
% (rc) Biomass half reaction : C(n)H(a)O(b)N(c) - Bacteria, Generic , N-Source : NH4+
% 0.208 CO2 + 0.0426 HCO3- + 0.0426 NH4+ + 1 H+ + 1 e- --> 0.2506 C'1'H'2.5'O'1'N'0.17 + 0.2932 H2O [ ?G = 27.4176 KJ/e-eq ]
%  
% Energy reaction : 
% + 0.0417 C6H12O6 + 0.25 H2O + 0.1667 N2 + 0.3333 H+ --> + 0.25 CO2 + 0.3333 NH4+  
%  
% Synthesis reaction : 
% + 0.0417 C6H12O6 + 0.0426 NH4+ + 0.0426 HCO3- --> 0.2506 C'1'H'2.5'O'1'N'0.17 + 0.042 CO2 + 0.0432 H2O  
%  
% Balanced equation using TEEM_2 : 
% [ fe = 0.87 ] [ fs = 0.13 ] [ e = 0.5 ]
% + 0.0417 C6H12O6 + 0.0055 HCO3- + 0.2123 H2O + 0.1452 N2 + 0.2904 H+ --> 0.0322 C'1'H'2.5'O'1'N'0.17 + 0.2232 CO2 + 0.285 NH4+  
%  
% Yield prediction :
% Yg/m = 25.469 [ grams_cells/mol_donor ]
% Yc/m = 0.774 [ mol_C_cells/mol_donor ]
% Yc/c = 0.129 [ mol_C_cell/mol_C_donor ]


a = [0.0417*6 0.0055  0.2123 0.1452 0.2904 0.0322 0.2232 0.285];  
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

% 
%  
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


% 
%  
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

%  
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



