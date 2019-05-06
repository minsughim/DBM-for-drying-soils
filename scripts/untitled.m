%analtic solution of pH while drying (change in concentration of net Z)

load('Standard_constanat_NH3_5ppb_HONO_1ppb.mat')
% at 25 degree
ListofZnet = logspace(-6,1,50);

for iZ = 1:50
        
        Znet = ListofZnet(iZ);        
        a1 = 1 + kHA*pNH3/Ka;
        a3 = -1*(Kw+KN*kHN*pHONO + K1c*kHC*pCO2);
        a4 = -2*K1c*K2C*kHC*pCO2;
        a2 = Znet;
        p = [a1 a2 a3 a4];
        r = roots(p);
        temp = max(r);
        pHsteady = -log10(temp);
        ResultListPH(iZ) = pHsteady;

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