load('Standard_constanat_NH3_5ppb_HONO_1ppb.mat')
% at 25 degree
NH3List = linspace(0,100,50);
HONOList = linspace(0,100,50);
for iNH3 = 1:50
    for iHONO = 1:50
        
        NH3ppb = NH3List(iNH3);
        HONOppb = HONOList(iHONO);
        
        pNH3temp = NH3ppb*pNH3/5;
        pHONOtemp = HONOppb*pHONO;
        
        a1 = 1 + kHA*pNH3temp/Ka;
        a3 = -1*(Kw+KN*kHN*pHONOtemp + K1c*kHC*pCO2);
        a4 = -2*K1c*K2C*kHC*pCO2;
        a2 = 0.1;
        p = [a1 a2 a3 a4];
        r = roots(p);
        temp = max(r);
        pHsteady = -log10(temp);
        ResultListNH3(iNH3,iHONO) = pNH3temp;
        ResultListHONO(iNH3,iHONO) = pHONOtemp;
        ResultListPH(iNH3,iHONO) = pHsteady;

        % a1H2ONO = 1 + kHA*pNH3/Ka +kHN*pHONO/KN2;
        % for i = 1:1000
        % a2 = a2List(i);
        % p = [a1H2ONO a2 a3 a4];
        % r = roots(p);
        % ListRootH2ONO(i,:) = r;
        % end
        % for i = 1:1000
        % temp = max(ListRootH2ONO(i,:));
        % pHsteadyH2ONO(i) = -log10(temp);
        % end
        %plot(a2List, pHsteady)
        %hold on
        %plot(a2List, pHsteadyH2ONO)
    end
end