load('Standard_constanat_NH3_5ppb_HONO_1ppb.mat')
% at 25 degree
%NH3List = logspace(-2,1,200);
%HONOList = logspace(-2,1,100);
pCO2 = 130*pCO2;
NH3List = 5;
HONOList = 1;
a = 0;
for iNH3 = 1:length(NH3List)
    for iHONO = 1:length(HONOList)
        
        NH3ppb = NH3List(iNH3);
        HONOppb = HONOList(iHONO);
        
        pNH3temp = NH3ppb*pNH3/5;
        pHONOtemp = HONOppb*pHONO;
        
        a1 = 1 + kHA*pNH3temp/Ka;
        a3 = -1*(Kw+KN*kHN*pHONOtemp + K1c*kHC*pCO2); 
        a4 = -2*K1c*K2C*kHC*pCO2;
        a2List = [-1*logspace(-9,1,100) logspace(-9,1,100)];
        a2List = sortrows(a2List');
        %a2List = [-0.1 -0.01 -0.001 0 0.001 0.01 0.1];
        for i = 1:length(a2List)
            a2 = a2List(i);
            p = [a1 a2 a3 a4];
            r = roots(p);
            ListRoot(i,:) = r;
            temp = max(r);
            pHsteady(i) = -log10(temp);
            a = a+1;
            ResultListNH3(iNH3,iHONO,i) = NH3ppb;
            ResultListHONO(iNH3,iHONO,i) = HONOppb;
            ResultListZnet(iNH3,iHONO,i) = a2List(i);
            ResultListPH(iNH3,iHONO,i) = pHsteady(i);
            
        end
        
        
        
        
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


temp = zeros(length(a2List),1);
temp2 = zeros(length(a2List),1);
temp(:) = ResultListZnet(1,1,:);
temp2(:) = ResultListPH(1,1,:);
plot(temp,temp2)
hold on