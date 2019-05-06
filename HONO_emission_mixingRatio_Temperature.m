numberOfSamples = 3;
EffluxHONO_honoMix_time = zeros(1728,numberOfSamples*3);
EffluxHONO_honoMix_WHC = zeros(288,numberOfSamples*3);
EffluxHONO_honoMix_WF = zeros(288,numberOfSamples*3);
EffluxNH3_honoMix_time = zeros(1728,numberOfSamples*3);
EffluxNH3_honoMix_WHC = zeros(288,numberOfSamples*3);

NO2Cons_time = zeros(1728,numberOfSamples*3);
NO2Prod_time = zeros(1728,numberOfSamples*3);
NO2tot_time = zeros(1728,numberOfSamples*3);
NO3Prod_time = zeros(1728,numberOfSamples*3);
NH4Cons_time = zeros(1728,numberOfSamples*3);

pHmax_time = zeros(1728,numberOfSamples*3);
pHmin_time = zeros(1728,numberOfSamples*3);
pHmean_time = zeros(1728,numberOfSamples*3);
pHstd_time = zeros(1728,numberOfSamples*3);

%% Varying NH3 and HONO at the same time
%ListHONOratio = [5 5 5 5; 0.1 5 10 15];
%ListHONOratio = [1 10 15 20; 1 1 1 1];
%ListHONOratio = [5 15 20 50;10 10 10 10];
%ListHONOratio = [5 5 5 5; 0.005 0.005 5 5;4 6 4 6];

%ListHONOratio = [20 20 20 20; 0.005 0.005 0.005 0.005;6 6 6 6; 25 27 28 30];
%ListHONOratio = [20 20 20 20 20 20; 0.005 0.005 0.005 0.005 0.005 0.005;6 6 6 6 6 6; 25 27 28 30 35 40];

%ListHONOratio = [5 5 5 5 5 5; 1 1 1 1 1 1;6 6 6 6 6 6; 25 27 28 30 32 34];
%ListHONOratio = [20 20 20 20 20 20; 0.005 0.005 0.005 0.005 0.005 0.005;6 6 6 6 6 6 ; 25 27 28 30 32 34];

ListHONOratio = [5 20 20 20 20 20;1 1 1 1 1 1;13 4 13 13 12 13;25 25 25 25 25 25];

%ListHONOratio = [20 20 20 20; 0.005 0.005 0.005 0.005;4 6 7 4; 30 30 30 15];


for iT= 1:numberOfSamples
% Various HONO mixing ratio
serial_id2 = sprintf('HT_BSC_photoY3_lineardecay_newT%d_NH3%d_HONO%d_Drying%d',ListHONOratio(4,iT), ListHONOratio(1,iT), ListHONOratio(2,iT), ListHONOratio(3,iT));
%serial_id2 = sprintf('HT_BSC_photoY3_Constant_dark_nobio_H2ONO_NH3%d_HONO%d_Drying%d',ListHONOratio(1,iT), ListHONOratio(2,iT), ListHONOratio(3,iT));
%serial_id2 = sprintf('HT_BSC_photoY3_lineardecay_NH3%d_HONO%d_Drying%d',ListHONOratio(1,iT), ListHONOratio(2,iT), ListHONOratio(3,iT));
load(strcat(serial_id2,'.mat'));
EffluxHONO_honoMix_time(:,3*(iT-1)+1) = timeLine -4*24;
EffluxNH3_honoMix_time(:,3*(iT-1)+1) = timeLine -4*24;
EffluxHONO_honoMix_WHC(:,3*(iT-1)+1) = avgWC(1153:1440)*100/0.24;
EffluxNH3_honoMix_WHC(:,3*(iT-1)+1) = avgWC(1153:1440)*100/0.24;
EffluxHONO_honoMix_WF(:,3*(iT-1)+1) = avgWF(1153:1440);


temp = mean(EffluxAll,3);
temp2 = std(EffluxAll,0,3);
EffluxHONO_honoMix_time(:,3*(iT-1)+2) = temp(:,4)*14/47;
EffluxHONO_honoMix_time(:,3*iT) = temp2(:,4)*14/47;
EffluxNH3_honoMix_time(:,3*(iT-1)+2) = temp(:,3)*14/17;
EffluxNH3_honoMix_time(:,3*iT) = temp2(:,3)*14/17;
EffluxHONO_honoMix_WHC(:,3*(iT-1)+2) = temp(1153:1440,4)*14/47;
EffluxHONO_honoMix_WHC(:,3*iT) = temp2(1153:1440,4)*14/47;
EffluxNH3_honoMix_WHC(:,3*(iT-1)+2) = temp(1153:1440,3)*14/17;
EffluxNH3_honoMix_WHC(:,3*iT) = temp2(1153:1440,3)*14/17;
EffluxHONO_honoMix_WF(:,3*(iT-1)+2) = temp(1153:1440,4)*14/47;
EffluxHONO_honoMix_WF(:,3*iT) = temp2(1153:1440,4)*14/47;


pHmax_time(:,3*(iT-1)+1) = timeLine -4*24;
pHmax_time(:,3*(iT-1)+2) = mean(maxpHtot');
pHmax_time(:,3*iT) = std(maxpHtot,0,2);
pHmin_time(:,3*(iT-1)+1) = timeLine -4*24;
pHmin_time(:,3*(iT-1)+2) = mean(minpHtot');
pHmin_time(:,3*iT) = std(minpHtot,0,2);
pHmean_time(:,3*(iT-1)+1) = timeLine -4*24;
pHmean_time(:,3*(iT-1)+2) = mean(meanpHtot');
pHmean_time(:,3*iT) = std(meanpHtot,0,2);
pHstd_time(:,3*(iT-1)+1) = timeLine -4*24;
pHstd_time(:,3*(iT-1)+2) = mean(stdpHtot');
pHstd_time(:,3*iT) = std(stdpHtot,0,2);


tempname = sprintf('stdpHtot_NH3%d_HONO%d_Drying%d', ListHONOratio(1,iT), ListHONOratio(2,iT), ListHONOratio(3,iT));
save(strcat(tempname,'.txt'),'stdpHtot','-ASCII')





% 
% figure(1)
% errorbar(EffluxHONO_honoMix_time(:,3*(iT-1)+1), EffluxHONO_honoMix_time(:,3*(iT-1)+2), EffluxHONO_honoMix_time(:,3*iT))
% hold on
% plot(EffluxHONO_honoMix_time(:,3*(iT-1)+1), EffluxHONO_honoMix_time(:,3*(iT-1)+2))
% errorbar(EffluxNH3_honoMix_time(:,3*(iT-1)+1), EffluxNH3_honoMix_time(:,3*(iT-1)+2), EffluxNH3_honoMix_time(:,3*iT))
% plot(EffluxNH3_honoMix_time(:,3*(iT-1)+1), EffluxNH3_honoMix_time(:,3*(iT-1)+2))
% figure(2)
% errorbar(EffluxHONO_honoMix_WHC(:,3*(iT-1)+1), EffluxHONO_honoMix_WHC(:,3*(iT-1)+2), EffluxHONO_honoMix_WHC(:,3*iT))
% plot(EffluxHONO_honoMix_WHC(:,3*(iT-1)+1), EffluxHONO_honoMix_WHC(:,3*(iT-1)+2))
% hold on
% plot(EffluxNH3_honoMix_WHC(:,3*(iT-1)+1), EffluxNH3_honoMix_WHC(:,3*(iT-1)+2))
% figure(3)
% errorbar(avgWF(1153:1440), temp(1153:1440,4)*14/47, temp2(1153:1440,4)*14/47)
% hold on
% plot(avgWF(1153:1440),temp(1153:1440,4)*14/47);
% figure(4)
% plot(EffluxNH3_honoMix_time(1153:1440,3*(iT-1)+1), EffluxHONO_honoMix_WHC(:,3*(iT-1)+1))
% hold on
% 
% temp = mean(NO2Cons,2);
% temp2 = std(NO2Cons,0,2);
% temp = NO2Cons(:,2);
% temp2 = zeros(1728,1);
% NO2Cons_time(:,3*(iT-1)+1) = timeLine -4*24;
% NO2Cons_time(:,3*(iT-1)+2) = temp*14/46;
% NO2Cons_time(:,3*iT) = temp2*14/46;
% 
% %temp = mean(NO2Prod,2);
% %temp2 = std(NO2Prod,0,2);
% temp = NO2Prod(:,2);
% temp2 = zeros(1728,1);
% NO2Prod_time(:,3*(iT-1)+1) = timeLine -4*24;
% NO2Prod_time(:,3*(iT-1)+2) = temp*14/46;
% NO2Prod_time(:,3*iT) = temp2*14/46;
% 
% %temp = mean(NO2tot,2);
% %temp2 = std(NO2tot,0,2);
% temp = NO2tot(:,2);
% temp2 = zeros(1728,1);
% NO2tot_time(:,3*(iT-1)+1) = timeLine -4*24;
% NO2tot_time(:,3*(iT-1)+2) = temp*14/46;
% NO2tot_time(:,3*iT) = temp2*14/46;
% 
% %temp = mean(NO3Prod,2);
% %temp2 = std(NO3Prod,0,2);
% temp = NO3Prod(:,2);
% temp2 = zeros(1728,1);
% NO3Prod_time(:,3*(iT-1)+1) = timeLine -4*24;
% NO3Prod_time(:,3*(iT-1)+2) = temp*14/62;
% NO3Prod_time(:,3*iT) = temp2*14/62;

% 
% %temp = mean(NH4Cons,2);
% %temp2 = std(NH4Cons,0,2);
% temp = NH4Cons(:,2);
% temp2 = zeros(1728,1);
% NH4Cons_time(:,3*(iT-1)+1) = timeLine -4*24;
% NH4Cons_time(:,3*(iT-1)+2) = temp*14/18;
% NH4Cons_time(:,3*iT) = temp2*14/18;

% 
% temp = zeros(300,4);
% for i =1:1728
% temp = AmmoniumProd(:,:,i).*totSoilWeight;
% AmmoniumConsSum(i) = sum(temp(:));
% temp = NitrateProd(:,:,i).*totSoilWeight;
% NitrateProdSum(i) = sum(temp(:));
% temp = NitriteProd(:,:,i);
% NitriteProdSum(i) = sum(temp(:));
% temp = NitriteCons(:,:,i);
% NitriteConsSum(i) = sum(temp(:));
% end
% 
% totalWeight = sum(totSoilWeight(:));
% 
% figure(5)
% plot(timeLine - 4*24, NitriteConsSum*10^9./totVerticalArea)
% hold on
% plot(timeLine - 4*24, NitriteProdSum*10^9./totVerticalArea)
% hold on
% plot(timeLine - 4*24, (NitriteProdSum+NitriteConsSum)*10^9./totVerticalArea)
% hold on
% 
% figure(6)
% plot(timeLine - 4*24, NitrateProdSum*10^9./totVerticalArea)
% hold on
% 
% figure(7)
% plot(timeLine - 4*24, AmmoniumConsSum*10^9./totVerticalArea)
% hold on
% 
% figure(8)
% plot(EffluxHONO_honoMix_WHC(:,3*(iT-1)+1),NitriteConsSum(1153:1440)*10^9./totVerticalArea);
% hold on
% plot(EffluxHONO_honoMix_WHC(:,3*(iT-1)+1),NitriteProdSum(1153:1440)*10^9./totVerticalArea);
% plot(EffluxHONO_honoMix_WHC(:,3*(iT-1)+1),(NitriteProdSum(1153:1440)+NitriteConsSum(1153:1440))*10^9./totVerticalArea);
% 


end
% 
% WaterContentsChange = zeros(288,4);
% WaterContentsChange(:,1) = timeLine(1153:1440) -4*24;
% WaterContentsChange(:,2) = EffluxNH3_honoMix_WHC(:,1);
% WaterContentsChange(:,3) = EffluxNH3_honoMix_WHC(:,4);
% WaterContentsChange(:,4) = EffluxNH3_honoMix_WHC(:,7);
% 
% save(strcat('./../../../../Fig/EffluxHONO_time_drying.txt'),'EffluxHONO_honoMix_time','-ASCII')
% save(strcat('./../../../../Fig/EffluxHONO_WHC_drying.txt'),'EffluxHONO_honoMix_WHC','-ASCII')
% save(strcat('./../../../../Fig/WaterContentsChange_drying.txt'),'WaterContentsChange','-ASCII')
% save(strcat('./../../../../Fig/NO2Cons_time_drying.txt'),'NO2Cons_time','-ASCII')
% save(strcat('./../../../../Fig/NO2Prod_time_drying.txt'),'NO2Prod_time','-ASCII')
% save(strcat('./../../../../Fig/NO2tot_time_drying.txt'),'NO2tot_time','-ASCII')
% save(strcat('./../../../../Fig/NO3Prod_time_drying.txt'),'NO3Prod_time','-ASCII')
% save(strcat('./../../../../Fig/NH4Cons_time_drying.txt'),'NH4Cons_time','-ASCII')


% cd ..
% load('Weber_2015_PNAS_Fig1.mat')
% figure(2)
% hold on
% plot(soilwatercontentWHC1,lightBSC1FNHONOfluxngm2s)
% plot(soilwatercontentWHC2,darkBSC32_2014FNHONOfluxngm2s)


% EffluxHONO_honoMix_time = zeros(1728,11);
% EffluxHONO_honoMix_WHC = zeros(288,11);
% WaterContentsChange = zeros(1728,2);
% 
% EffluxNH3_honoMix_time = zeros(1728,11);
% EffluxNH3_honoMix_WHC = zeros(288,11);
% 
% ListHONOratio = [0.1 1 5 10 15];
% 
% for iT= 1:5
% % Various HONO mixing ratio
% serial_id2 = sprintf('HT_BSC_photoY3_lineardecay_HONO%d_Drying%d', ListHONOratio(iT), 4);
% load(strcat(serial_id2,'.mat'));
% EffluxHONO_honoMix_time(:,1) = timeLine -4*24;
% EffluxNH3_honoMix_time(:,1) = timeLine -4*24;
% EffluxHONO_honoMix_WHC(:,1) = avgWC(1153:1440)*100/0.28;
% EffluxNH3_honoMix_WHC(:,1) = avgWC(1153:1440)*100/0.28;
% 
% temp = mean(EffluxAll,3);
% temp2 = std(EffluxAll,0,3);
% EffluxHONO_honoMix_time(:,2*iT) = temp(:,4)*14/47;
% EffluxHONO_honoMix_time(:,2*iT+1) = temp2(:,4)*14/47;
% EffluxNH3_honoMix_time(:,2*iT) = temp(:,3)*14/17;
% EffluxNH3_honoMix_time(:,2*iT+1) = temp2(:,3)*14/17;
% EffluxHONO_honoMix_WHC(:,2*iT) = temp(1153:1440,4)*14/47;
% EffluxHONO_honoMix_WHC(:,2*iT+1) = temp2(1153:1440,4)*14/47;
% EffluxNH3_honoMix_WHC(:,2*iT) = temp(1153:1440,3)*14/17;
% EffluxNH3_honoMix_WHC(:,2*iT+1) = temp2(1153:1440,3)*14/17;
% 
% figure(1)
% errorbar(EffluxHONO_honoMix_time(:,1), EffluxHONO_honoMix_time(:,2*iT), EffluxHONO_honoMix_time(:,2*iT+1))
% hold on
% plot(EffluxHONO_honoMix_time(:,1), EffluxHONO_honoMix_time(:,2*iT))
% errorbar(EffluxNH3_honoMix_time(:,1), EffluxNH3_honoMix_time(:,2*iT), EffluxNH3_honoMix_time(:,2*iT+1))
% plot(EffluxNH3_honoMix_time(:,1), EffluxNH3_honoMix_time(:,2*iT))
% figure(2)
% plot(EffluxHONO_honoMix_WHC(:,1), EffluxHONO_honoMix_WHC(:,2*iT))
% hold on
% plot(EffluxNH3_honoMix_WHC(:,1), EffluxNH3_honoMix_WHC(:,2*iT))
% 
% end
% 
% WaterContentsChange(:,1) = timeLine -4*24;
% WaterContentsChange(:,2) = avgWC;
% 
% save(strcat('./../../../../Fig/WaterContentsChange.txt'),'WaterContentsChange','-ASCII')
% save(strcat('./../../../../Fig/EffluxHONO_honoMx_time.txt'),'EffluxHONO_honoMix_time','-ASCII')
% save(strcat('./../../../../Fig/EffluxNH3_honoMix_time.txt'),'EffluxNH3_honoMix_time','-ASCII')
% save(strcat('./../../../../Fig/EffluxHONO_honoMix_WHC.txt'),'EffluxHONO_honoMix_WHC','-ASCII')
% save(strcat('./../../../../Fig/EffluxNH3_honoMix_WHC.txt'),'EffluxNH3_honoMix_WHC','-ASCII')

save(strcat('EffluxHONO_ext_time2.txt'),'EffluxHONO_honoMix_time','-ASCII')
save(strcat('EffluxNH3_ext_time2.txt'),'EffluxNH3_honoMix_time','-ASCII')
save(strcat('EffluxHONO_ext_WHC2.txt'),'EffluxHONO_honoMix_WHC','-ASCII')
save(strcat('EffluxNH3_ext_WHC2.txt'),'EffluxNH3_honoMix_WHC','-ASCII')
% 
save(strcat('pHmax_time.txt'),'pHmax_time','-ASCII')
save(strcat('pHmin_time.txt'),'pHmin_time','-ASCII')
save(strcat('pHmean_time.txt'),'pHmean_time','-ASCII')
save(strcat('pHstd_time.txt'),'pHstd_time','-ASCII')
%

% EffluxHONO_nh3Mix_time = zeros(1728,11);
% EffluxHONO_nh3Mix_WHC = zeros(288,11);
% EffluxNH3_nh3Mix_time = zeros(1728,11);
% EffluxNH3_nh3Mix_WHC = zeros(288,11);
% 
% ListNH3ratio = [1 5 10 15 20];
% 
% for iT= 1:5
% % Various HONO mixing ratio
% serial_id2 = sprintf('HT_BSC_photoY3_lineardecay_NH3%d_Drying%d', ListNH3ratio(iT), 4);
% load(strcat(serial_id2,'.mat'));
% EffluxHONO_nh3Mix_time(:,1) = timeLine -4*24;
% EffluxNH3_nh3Mix_time(:,1) = timeLine -4*24;
% EffluxHONO_nh3Mix_WHC(:,1) = avgWC(1153:1440)*100/0.28;
% EffluxNH3_nh3Mix_WHC(:,1) = avgWC(1153:1440)*100/0.28;
% 
% temp = mean(EffluxAll,3);
% temp2 = std(EffluxAll,0,3);
% EffluxHONO_nh3Mix_time(:,2*iT) = temp(:,4)*14/47;
% EffluxHONO_nh3Mix_time(:,2*iT+1) = temp2(:,4)*14/47;
% EffluxNH3_nh3Mix_time(:,2*iT) = temp(:,3)*14/17;
% EffluxNH3_nh3Mix_time(:,2*iT+1) = temp2(:,3)*14/17;
% EffluxHONO_nh3Mix_WHC(:,2*iT) = temp(1153:1440,4)*14/47;
% EffluxHONO_nh3Mix_WHC(:,2*iT+1) = temp2(1153:1440,4)*14/47;
% EffluxNH3_nh3Mix_WHC(:,2*iT) = temp(1153:1440,3)*14/17;
% EffluxNH3_nh3Mix_WHC(:,2*iT+1) = temp2(1153:1440,3)*14/17;
% 
% figure(1)
% errorbar(EffluxHONO_nh3Mix_time(:,1), EffluxHONO_nh3Mix_time(:,2*iT), EffluxHONO_nh3Mix_time(:,2*iT+1))
% hold on
% plot(EffluxHONO_nh3Mix_time(:,1), EffluxHONO_nh3Mix_time(:,2*iT))
% errorbar(EffluxNH3_nh3Mix_time(:,1), EffluxNH3_nh3Mix_time(:,2*iT), EffluxNH3_nh3Mix_time(:,2*iT+1))
% plot(EffluxNH3_nh3Mix_time(:,1), EffluxNH3_nh3Mix_time(:,2*iT))
% figure(2)
% plot(EffluxHONO_nh3Mix_WHC(:,1), EffluxHONO_nh3Mix_WHC(:,2*iT))
% hold on
% plot(EffluxNH3_nh3Mix_WHC(:,1), EffluxNH3_nh3Mix_WHC(:,2*iT))
% 
% end
% 
% save(strcat('./../../../../Fig/EffluxHONO_nh3Mix_time.txt'),'EffluxHONO_nh3Mix_time','-ASCII')
% save(strcat('./../../../../Fig/EffluxNH3_nh3Mix_time.txt'),'EffluxNH3_nh3Mix_time','-ASCII')
% save(strcat('./../../../../Fig/EffluxHONO_nh3Mix_WHC.txt'),'EffluxHONO_nh3Mix_WHC','-ASCII')
% save(strcat('./../../../../Fig/EffluxNH3_nh3Mix_WHC.txt'),'EffluxNH3_nh3Mix_WHC','-ASCII')
% 
% 


