function [sitesC2t,sitesC3, sitesC5, sitesC6,sitesC7, sitesC22,sitesC23, sitesC24,sitesC25 ,sitesC26, sitesC27, sitesCaTemp] = PHestimationMex(sitesC2t,sitesC3, sitesC5, sitesC6, sitesC7, sitesC22,sitesC23, sitesC24,sitesC25,sitesC26,sitesC27,sitesCaTemp, DeltaT, ZionTemp,pKTotList,thetaTC)
%N1 = O2
%Carbon related
%N2 = CO2
%N3 = HCO3-
%N4 = carbohydrates (CH2O) : sugar like glucose
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
% 9: CH1.8O0.5N0.2                 32.91
%10: HCO3                          61.02
%11: NH4                           18.04
%12: CO3                           60.01
%13: NO2-                          46.01
%14: HONO                          47.013
Molwt = [30.0263 31.9988 17.0306 62.0049 28.0135 24.6263 44.0098 18.0153 32.9114 61.0171 18.0385 60.0092 46.01 47.013];
%change all units based on kmol/m^3

[m,n] = size(sitesC6);
sitesTemp = zeros(m,n,11);

sitesTemp(:,:,1) = sitesC6*10^(-3)/Molwt(11);%kmol/m^3 (M) change all unites
%sitesTemp(:,:,1) = (sitesC6-sitesC23)*10^(-3)/Molwt(11);%kmol/m^3 (M) change all unites
sitesTemp(:,:,2) = sitesC23*10^(-3)/Molwt(3);
sitesTemp(:,:,3) = sitesC2t*10^(-3)/Molwt(7);
sitesTemp(:,:,4) = sitesC3*10^(-3)/Molwt(10);
sitesTemp(:,:,5) = sitesC22*10^(-3)/Molwt(12);
sitesTemp(:,:,6) = sitesC7*10^(-3)/Molwt(13); %nitrite
sitesTemp(:,:,7) = sitesC26*10^(-3)/Molwt(14); %HONO
sitesTemp(:,:,8) = sitesCaTemp(:,:,1); %CaCO3 complex in M
sitesTemp(:,:,9) = sitesCaTemp(:,:,2); %CaCO3 precipitation in M
sitesTemp(:,:,10) = sitesC27; %Ca in M
sitesTemp(:,:,11) = sitesC24; %H in M

ZTemp = ZionTemp -sitesC5*10^(-3)/Molwt(4); %nitrate came infor charge balance for the mother material of the soil, here other cation and anion can be added.

%M = [eye(5) zeros(5,1); zeros(1,6)];
%options = odeset('AbsTol',1E-6,'RelTol',1E-6,'Mass',M, 'event', @eventPH);
options = odeset('AbsTol',1E-8,'RelTol',1E-8,'NonNegative',1:10);%, 'event', @eventPH);
%options = odeset('AbsTol',1E-10,'RelTol',1E-6);

y0 = zeros(10,1);
pKs = zeros(6,1);
ions = [1 0 0 -1 -2 -1 0 0 0 2];

for i = 1:m
    for j = 1:n
        pKs(:) = pKTotList(i,j,:);
        y0(:) = sitesTemp(i,j,1:10);
        y0(y0<=0) = 1.E-16; % to inhibit the divergence in the code, small number was assgined for zero concentrations (caused by fast reaction term)
        try            
            [T2,Y2] = pH_mex_C_3(DeltaT, y0, [ZTemp(i,j) thetaTC(i,j) pKs'], [],[1e-8 1e-8 10]);
            solY = Y2(:,end);
        catch
            %warning(lasterr);
            [T,Y] = ode15s(@(t,y2) HONO_forced_kinetics_Tdep_calcite(y2,ZTemp(i,j),pKs,thetaTC(i,j)), [0 DeltaT], y0,options);
            solY = Y(end,:)';
        end    
        
        solY(solY<=0) = 1.E-16; % force non-negative solution
        totIon = ions*solY +ZTemp(i,j); %Other cation and anoion included here for charge balancing
        y11 = 0.5*(sqrt(totIon*totIon+4*10^(-14))-totIon);
        
        sitesC6(i,j) = solY(1)*10^(3)*Molwt(11);%kmol/m^3 (M) change all unites
        sitesC23(i,j) = solY(2)*10^(3)*Molwt(3);
        sitesC2t(i,j) = solY(3)*10^(3)*Molwt(7);
        sitesC3(i,j) = solY(4)*10^(3)*Molwt(10);
        sitesC22(i,j) = solY(5)*10^(3)*Molwt(12);
        sitesC24(i,j) = y11;
        sitesC25(i,j) = -log10(y11);
        sitesC7(i,j) = solY(6)*10^(3)*Molwt(13);
        sitesC26(i,j) = solY(7)*10^(3)*Molwt(14);
        sitesCaTemp(i,j,1) = solY(8); %CaCO3 complex in M
        sitesCaTemp(i,j,2) = solY(9); %CaCO3 precipitation in M
        sitesC27(i,j) = solY(10); %Ca in M
        
    end
end



