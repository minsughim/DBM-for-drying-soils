
numberOFP = 8;
g=graph('2d-deg6-trilattice',{n,m,Lp});
listV = g.vertices;
listE = g.edges;
%periodic boundary condition for y direction
listE(m*n-1,m) = 1;
listE(m*n,m) = 1;
listE(m,m*n-1) = 1;
listE(m,m*n) = 1;
for i = 1:2:(m-1)
    listE((n-1)*m+i ,i) = 1;
    listE(i,(n-1)*m+i) = 1;
end
for i = 2:2:(m-2)
    listE((n-1)*m+i-1 ,i) = 1;
    listE((n-1)*m+i ,i) = 1;
    listE((n-1)*m+i+1 ,i) = 1;
    listE(i,(n-1)*m+i-1) = 1;
    listE(i,(n-1)*m+i) = 1;
    listE(i,(n-1)*m+i+1) = 1;
end

positionListT = zeros(m,n,2);
positionListT(:,:,1) = reshape(listV(:,1), [m,n])*1000;
positionListT(:,:,2) = reshape(listV(:,2), [m,n])*(-1000);
vmin = 1;
vmax = numberOFP;
dv = vmax-vmin;
colormatirx = ones(3,numberOFP);
colormatirx(:,1) = [0.47,0.67,0.19];%green : phototrophs 
colormatirx(:,2) = [0.47,0.67,0.19];%green : phototrophs 
colormatirx(:,3) = [0.47,0.67,0.19];%green : phototrophs 
colormatirx(:,4) = [0.47,0.67,0.19];%green : phototrophs 
colormatirx(:,5) = [0.93,0.69,0.13];%Yellow : aerobes
colormatirx(:,6) = [0.49,0.18,0.56];%purple : anaerobes
colormatirx(:,7) = [0,0.45,0.74];%blue : nitrifier
colormatirx(:,8) = [0,0.45,0.74];%blue : nitrifier




%% For chemical profile without cells
%
% HoriList = 0:Lp:Lp*(n-1);
% figure(2)
% subplot(1,2,1)
% contour(xlist,ylist,  timeConcDist{1}(:,:,dayT))
% caxis([0 35])
% colorbar;
% subplot(1,2,2)
% contour(xlist,ylist,  timeConcDist{1}(:,:,nightT))
% caxis([0 35])
% colorbar;
%
% xlist = reshape(listV(:,1), [m,n])*1000;
% ylist = reshape(listV(:,2), [m,n])*(-1000);
% LocalGasContentstext = zeros(2000,3);
% LocalGasContentstext(:,1) = xlist(:);
% LocalGasContentstext(:,2) = ylist(:);
% LocalGasContentstext(:,3) = LocalGasContents(:);
% LocalGasContentstext2 = sortrows(LocalGasContentstext,-2);
% fileID = fopen('./../../../Fig/LocalGasContents.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',LocalGasContentstext2(i,1),LocalGasContentstext2(i,2)*(-1),LocalGasContentstext2(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
% O2daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = sitesC{1}(:,:);
% O2daytext(:,1) = xlist(:);
% O2daytext(:,2) = ylist(:);
% O2daytext(:,3) = temp(:);
% O2daytext = sortrows(O2daytext,-2);
% fileID = fopen('./../../../Fig/O2profile.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% pHdaytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = sitesC2{5}(:,:);
% pHdaytext(:,1) = xlist(:);
% pHdaytext(:,2) = ylist(:);
% pHdaytext(:,3) = temp(:);
% pHdaytext = sortrows(pHdaytext,-2);
% fileID = fopen('./../../../Fig/pHprofile.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',pHdaytext(i,1),pHdaytext(i,2)*(-1),pHdaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
% TANdaytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = sitesC{6}(:,:)+sitesC2{3}(:,:);
% TANdaytext(:,1) = xlist(:);
% TANdaytext(:,2) = ylist(:);
% TANdaytext(:,3) = temp(:);
% TANdaytext = sortrows(TANdaytext,-2);
% fileID = fopen('./../../../Fig/TANprofile.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',TANdaytext(i,1),TANdaytext(i,2)*(-1),TANdaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% NO3daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = sitesC{5}(:,:);
% NO3daytext(:,1) = xlist(:);
% NO3daytext(:,2) = ylist(:);
% NO3daytext(:,3) = temp(:);
% NO3daytext = sortrows(NO3daytext,-2);
% fileID = fopen('./../../../Fig/NO3profile.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% Islanddaytext = zeros(2000,3);
% %temp = zeros(100,20);
% %temp = InvasedIsland;
% Islanddaytext(:,1) = xlist(:);
% Islanddaytext(:,2) = ylist(:);
% Islanddaytext(:,3) = temp(:);
% Islanddaytext = sortrows(Islanddaytext,-2);
% fileID = fopen('./../../../Fig/Islandprofile.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',Islanddaytext(i,1),Islanddaytext(i,2)*(-1),Islanddaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%


% %% For chemical profile with cells at water filled capacity
% %
% nightT = 1368;
% dayT = 1224;
% HoriList = 0:Lp:Lp*(n-1);
% % figure(2)
% % subplot(1,2,1)
% % contour(xlist,ylist,  timeConcDist{1}(:,:,dayT))
% % caxis([0 35])
% % colorbar;
% % subplot(1,2,2)
% % contour(xlist,ylist,  timeConcDist{1}(:,:,nightT))
% % caxis([0 35])
% % colorbar;
%
% O2daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{1}(:,:,dayT);
% O2daytext(:,1) = xlist(:);
% O2daytext(:,2) = ylist(:);
% O2daytext(:,3) = temp(:);
% O2daytext = sortrows(O2daytext,-2);
% fileID = fopen('./../../../Fig/O2profileDay.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% O2daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{1}(:,:,nightT);
% O2daytext(:,1) = xlist(:);
% O2daytext(:,2) = ylist(:);
% O2daytext(:,3) = temp(:);
% O2daytext = sortrows(O2daytext,-2);
% fileID = fopen('./../../../Fig/O2profileNight.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
%
% pHdaytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist2{5}(:,:,dayT);
% pHdaytext(:,1) = xlist(:);
% pHdaytext(:,2) = ylist(:);
% pHdaytext(:,3) = temp(:);
% pHdaytext = sortrows(pHdaytext,-2);
% fileID = fopen('./../../../Fig/pHprofileDay.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',pHdaytext(i,1),pHdaytext(i,2)*(-1),pHdaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
% pHdaytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist2{5}(:,:,nightT);
% pHdaytext(:,1) = xlist(:);
% pHdaytext(:,2) = ylist(:);
% pHdaytext(:,3) = temp(:);
% pHdaytext = sortrows(pHdaytext,-2);
% fileID = fopen('./../../../Fig/pHprofileNight.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',pHdaytext(i,1),pHdaytext(i,2)*(-1),pHdaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% TANdaytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{6}(:,:,dayT)+timeConcDist2{3}(:,:,dayT);
% TANdaytext(:,1) = xlist(:);
% TANdaytext(:,2) = ylist(:);
% TANdaytext(:,3) = temp(:);
% TANdaytext = sortrows(TANdaytext,-2);
% fileID = fopen('./../../../Fig/TANprofileDay.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',TANdaytext(i,1),TANdaytext(i,2)*(-1),TANdaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
% TANdaytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{6}(:,:,nightT)+timeConcDist2{3}(:,:,nightT);
% TANdaytext(:,1) = xlist(:);
% TANdaytext(:,2) = ylist(:);
% TANdaytext(:,3) = temp(:);
% TANdaytext = sortrows(TANdaytext,-2);
% fileID = fopen('./../../../Fig/TANprofileNight.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',TANdaytext(i,1),TANdaytext(i,2)*(-1),TANdaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
%
% NO3daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{5}(:,:,dayT);
% NO3daytext(:,1) = xlist(:);
% NO3daytext(:,2) = ylist(:);
% NO3daytext(:,3) = temp(:);
% NO3daytext = sortrows(NO3daytext,-2);
% fileID = fopen('./../../../Fig/NO3profileDay.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% NO3daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{5}(:,:,nightT);
% NO3daytext(:,1) = xlist(:);
% NO3daytext(:,2) = ylist(:);
% NO3daytext(:,3) = temp(:);
% NO3daytext = sortrows(NO3daytext,-2);
% fileID = fopen('./../../../Fig/NO3profileNight.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
%
% NO3daytext = zeros(2000,3);
% temp = zeros(100,20);
% temp = timeConcDist{5}(:,:,nightT);
% NO3daytext(:,1) = xlist(:);
% NO3daytext(:,2) = ylist(:);
% NO3daytext(:,3) = temp(:);
% NO3daytext = sortrows(NO3daytext,-2);
% fileID = fopen('./../../../Fig/NO3profileNight.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%


%% Microbial organism distribution
Rad3Over2 = sqrt(3)/2;
Rad3 = sqrt(3);

Lp = (depthList(2)-depthList(1))/Rad3Over2; % size of patch for toplayer [m]
dl = Lp/Rad3; % length of the edge of hexagon (top layer)
patchArea = Rad3Over2*dl*dl;
totSoilVolume = patchArea*TotalporeD./porosityM;
totSoilWeight = totSoilVolume*1300000;


timeT = 1368;

for iSp = 1:8;
    temp1 = full(popMovie{iSp,timeT});
    temp3 = BioMassDist{iSp,timeT};
    temp2 = full(ActivityCells{iSp,timeT});
    temp = temp2*popWalkers(timeT,iSp)./temp1.*temp3.*10^6*3600./totSoilWeight;
    temp(isnan(temp)) = 0;
    activityProf(:,:,iSp) = temp(:,:);
end

%% Plotting for Python
serial_id = sprintf('bioCrust_day_pop');
popTemp = activityProf*10^7;
tempListPop1 = 0;
for i=1:(m*0.5)
    for j = 1:n
        x = positionListT(i,j,1);
        y = positionListT(i,j,2);
        
        for ik = 1:8
            colorV = colormatirx(:,ik);
            r2 =log10(popTemp(i,j,ik));                
            r1 = 1*(r2>1);
            if r1 ~= 0
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,1) = x; %in cm scale
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,2) = y;                               
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,3) = colorV(1);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,4) = colorV(2);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,5) = colorV(3);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,6) = (r2>7)*1 +(r2/14)*(r2<=7);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,7) = i; %in cm scale
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,8) = j;
                tempListPop1 = tempListPop1 + r1;
            end            
        end
    end
end

totalPop = length(ListOfWalkersText1);
list = randperm(totalPop);
plotListWalkers = zeros(totalPop,8);
for i = 1:totalPop
    tempL = list(i);
    plotListWalkers(i,:) = ListOfWalkersText1(tempL,:);
end

save(strcat('./../../../Python/',serial_id,'.txt'),'plotListWalkers','-ASCII')

clear ListOfWalkersText1



PHOTOdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = activityProf(:,:,1)+activityProf(:,:,2)+activityProf(:,:,3)+activityProf(:,:,4);
PHOTOdaytext(:,1) = xlist(:);
PHOTOdaytext(:,2) = ylist(:);
PHOTOdaytext(:,3) = temp(:);
PHOTOdaytext = sortrows(PHOTOdaytext,-2);
fileID = fopen('./../../../Fig/PHOTOprofileNight.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',PHOTOdaytext(i,1),PHOTOdaytext(i,2)*(-1),PHOTOdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/PHOTO_SpatialAvNight.txt'),'SpatialAv','-ASCII')


AEROdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,5);
AEROdaytext(:,1) = xlist(:);
AEROdaytext(:,2) = ylist(:);
AEROdaytext(:,3) = temp(:);
AEROdaytext = sortrows(AEROdaytext,-2);
fileID = fopen('./../../../Fig/AEROprofileNight.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',AEROdaytext(i,1),AEROdaytext(i,2)*(-1),AEROdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/AERO_SpatialAvNight.txt'),'SpatialAv','-ASCII')

ANAEROdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,6);
ANAEROdaytext(:,1) = xlist(:);
ANAEROdaytext(:,2) = ylist(:);
ANAEROdaytext(:,3) = temp(:);
ANAEROdaytext = sortrows(ANAEROdaytext,-2);
fileID = fopen('./../../../Fig/ANAEROprofileNight.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',ANAEROdaytext(i,1),ANAEROdaytext(i,2)*(-1),ANAEROdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/ANAERO_SpatialAvNight.txt'),'SpatialAv','-ASCII')


AOBdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,7);
AOBdaytext(:,1) = xlist(:);
AOBdaytext(:,2) = ylist(:);
AOBdaytext(:,3) = temp(:);
AOBdaytext = sortrows(AOBdaytext,-2);
fileID = fopen('./../../../Fig/AOBprofileNight.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',AOBdaytext(i,1),AOBdaytext(i,2)*(-1),AOBdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/AOB_SpatialAvNight.txt'),'SpatialAv','-ASCII')


NOBdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,8);
NOBdaytext(:,1) = xlist(:);
NOBdaytext(:,2) = ylist(:);
NOBdaytext(:,3) = temp(:);
NOBdaytext = sortrows(NOBdaytext,-2);
fileID = fopen('./../../../Fig/NOBprofileNight.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',NOBdaytext(i,1),NOBdaytext(i,2)*(-1),NOBdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/NOB_SpatialAvNight.txt'),'SpatialAv','-ASCII')

Nitridaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,7)+activityProf(:,:,8);
Nitridaytext(:,1) = xlist(:);
Nitridaytext(:,2) = ylist(:);
Nitridaytext(:,3) = temp(:);
Nitridaytext = sortrows(Nitridaytext,-2);
fileID = fopen('./../../../Fig/NitriprofileNight.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',Nitridaytext(i,1),Nitridaytext(i,2)*(-1),Nitridaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/Nitri_SpatialAvNight.txt'),'SpatialAv','-ASCII')


timeT = 1224;
for iSp = 1:8;
    temp1 = full(popMovie{iSp,timeT});
    temp3 = BioMassDist{iSp,timeT};
    temp2 = full(ActivityCells{iSp,timeT});
    temp = temp2*popWalkers(timeT,iSp)./temp1.*temp3.*10^6*3600./totSoilWeight;
    temp(isnan(temp)) = 0;
    activityProf(:,:,iSp) = temp(:,:);
end

serial_id = sprintf('bioCrust_night_pop');
popTemp = activityProf*10^7;
tempListPop1 = 0;
for i=1:(m*0.5)
    for j = 1:n
        x = positionListT(i,j,1);
        y = positionListT(i,j,2);
        
        for ik = 1:8
            colorV = colormatirx(:,ik);
            r2 =log10(popTemp(i,j,ik));                
            r1 = 1*(r2>1);
            if r1 ~= 0
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,1) = x; %in cm scale
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,2) = y;                               
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,3) = colorV(1);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,4) = colorV(2);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,5) = colorV(3);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,6) = (r2>7)*1 +(r2/14)*(r2<=7);
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,7) = i; %in cm scale
                ListOfWalkersText1(tempListPop1+1:tempListPop1 +r1,8) = j;
                tempListPop1 = tempListPop1 + r1;
            end            
        end
    end
end

totalPop = length(ListOfWalkersText1);
list = randperm(totalPop);
plotListWalkers = zeros(totalPop,8);
for i = 1:totalPop
    tempL = list(i);
    plotListWalkers(i,:) = ListOfWalkersText1(tempL,:);
end
save(strcat('./../../../Python/',serial_id,'.txt'),'plotListWalkers','-ASCII')

clear ListOfWalkersText1



PHOTOdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = activityProf(:,:,1)+activityProf(:,:,2)+activityProf(:,:,3)+activityProf(:,:,4);
PHOTOdaytext(:,1) = xlist(:);
PHOTOdaytext(:,2) = ylist(:);
PHOTOdaytext(:,3) = temp(:);
PHOTOdaytext = sortrows(PHOTOdaytext,-2);

fileID = fopen('./../../../Fig/PHOTOprofileDay.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',PHOTOdaytext(i,1),PHOTOdaytext(i,2)*(-1),PHOTOdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/PHOTO_SpatialAvDay.txt'),'SpatialAv','-ASCII')


AEROdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,5);
AEROdaytext(:,1) = xlist(:);
AEROdaytext(:,2) = ylist(:);
AEROdaytext(:,3) = temp(:);
AEROdaytext = sortrows(AEROdaytext,-2);
fileID = fopen('./../../../Fig/AEROprofileDay.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',AEROdaytext(i,1),AEROdaytext(i,2)*(-1),AEROdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/AERO_SpatialAvDay.txt'),'SpatialAv','-ASCII')

ANAEROdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,6);
ANAEROdaytext(:,1) = xlist(:);
ANAEROdaytext(:,2) = ylist(:);
ANAEROdaytext(:,3) = temp(:);
ANAEROdaytext = sortrows(ANAEROdaytext,-2);
fileID = fopen('./../../../Fig/ANAEROprofileDay.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',ANAEROdaytext(i,1),ANAEROdaytext(i,2)*(-1),ANAEROdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/ANAERO_SpatialAvDay.txt'),'SpatialAv','-ASCII')


AOBdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,7);
AOBdaytext(:,1) = xlist(:);
AOBdaytext(:,2) = ylist(:);
AOBdaytext(:,3) = temp(:);
AOBdaytext = sortrows(AOBdaytext,-2);
fileID = fopen('./../../../Fig/AOBprofileDay.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',AOBdaytext(i,1),AOBdaytext(i,2)*(-1),AOBdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/AOB_SpatialAvDay.txt'),'SpatialAv','-ASCII')


NOBdaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,8);
NOBdaytext(:,1) = xlist(:);
NOBdaytext(:,2) = ylist(:);
NOBdaytext(:,3) = temp(:);
NOBdaytext = sortrows(NOBdaytext,-2);
fileID = fopen('./../../../Fig/NOBprofileDay.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',NOBdaytext(i,1),NOBdaytext(i,2)*(-1),NOBdaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/NOB_SpatialAvDay.txt'),'SpatialAv','-ASCII')

Nitridaytext = zeros(2000,3);
temp = zeros(100,20);
temp = activityProf(:,:,7)+activityProf(:,:,8);
Nitridaytext(:,1) = xlist(:);
Nitridaytext(:,2) = ylist(:);
Nitridaytext(:,3) = temp(:);
Nitridaytext = sortrows(Nitridaytext,-2);
fileID = fopen('./../../../Fig/NitriprofileDay.txt','w');
for i = 1:2000
    fprintf(fileID,'%f %f %f\n',Nitridaytext(i,1),Nitridaytext(i,2)*(-1),Nitridaytext(i,3));
    if rem(i,20) == 0
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
SpatialAv = zeros(100,3);
SpatialAv(:,1) = depthList;
SpatialAv(:,2) = mean(temp,2);
SpatialAv(:,3)  = std(temp,0,2);
save(strcat('./../../../Fig/Nitri_SpatialAvDay.txt'),'SpatialAv','-ASCII')

%
% ActivityCells = cell(5,1);
%
% for timeT = 1:1339
%     ActivityCells{1}(:,:,timeT) = full(Mumean(timeT,1)*popMovie{1,timeT}+Mumean(timeT,2)*popMovie{2,timeT}+Mumean(timeT,3)*popMovie{3,timeT}+Mumean(timeT,4)*popMovie{4,timeT})*rho*Vu*10^9*3600./totSoilWeight;
%     ActivityCells{2}(:,:,timeT) = full(Mumean(timeT,5)*popMovie{5,timeT})*rho*Vu*10^9*3600./totSoilWeight;
%     ActivityCells{3}(:,:,timeT) = full(Mumean(timeT,6)*popMovie{6,timeT})*rho*Vu*10^9*3600./totSoilWeight;
%     ActivityCells{4}(:,:,timeT) = full(Mumean(timeT,7)*popMovie{7,timeT})*rho*Vu*10^9*3600./totSoilWeight;
%     ActivityCells{5}(:,:,timeT) = full(Mumean(timeT,8)*popMovie{8,timeT})*rho*Vu*10^9*3600./totSoilWeight;
% end


%
% Islanddaytext = zeros(2000,3);
% %temp = zeros(100,20);
% %temp = InvasedIsland;
% Islanddaytext(:,1) = xlist(:);
% Islanddaytext(:,2) = ylist(:);
% Islanddaytext(:,3) = temp(:);
% Islanddaytext = sortrows(Islanddaytext,-2);
% fileID = fopen('./../../../Fig/Islandprofile.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',Islanddaytext(i,1),Islanddaytext(i,2)*(-1),Islanddaytext(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%




%
% LocalGasContentstext = zeros(2000,3);
% LocalGasContentstext(:,1) = xlist(:);
% LocalGasContentstext(:,2) = ylist(:);
% LocalGasContentstext(:,3) = LocalGasContents(:);
% LocalGasContentstext2 = sortrows(LocalGasContentstext,-2);
% fileID = fopen('./../../../Fig/LocalGasContents.txt','w');
% for i = 1:2000
% fprintf(fileID,'%f %f %f\n',LocalGasContentstext2(i,1),LocalGasContentstext2(i,2)*(-1),LocalGasContentstext2(i,3));
% if rem(i,20) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);




%
% waterfilmtext = zeros(800,3);
% waterfilmtext(:,1) = xlist(:);
% waterfilmtext(:,2) = ylist(:);
% waterfilmtext(:,3) = waterfilm(:);
% waterfilmtext2 = sortrows(waterfilmtext,-2);
% fileID = fopen('./../../../Fig/waterfilm.txt','w');
% for i = 1:800
% fprintf(fileID,'%f %f %f\n',waterfilmtext2(i,1),waterfilmtext2(i,2)*(-1),waterfilmtext2(i,3)*1000000);
% if rem(i,8) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
% O2daytext = zeros(800,3);
% temp = zeros(100,8);
% temp = timeConcDist{1}(:,:,dayT);
% O2daytext(:,1) = xlist(:);
% O2daytext(:,2) = ylist(:);
% O2daytext(:,3) = temp(:);
%
% O2daytext = sortrows(O2daytext,-2);
% fileID = fopen('./../../../Fig/O2day.txt','w');
% for i = 1:800
% fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
% if rem(i,8) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);
%
% O2nighttext = zeros(800,3);
% temp = zeros(100,8);
% temp = timeConcDist{1}(:,:,nightT);
% O2nighttext(:,1) = xlist(:);
% O2nighttext(:,2) = ylist(:);
% O2nighttext(:,3) = temp(:);
%
% O2nighttext = sortrows(O2nighttext,-2);
% fileID = fopen('./../../../Fig/O2night.txt','w');
% for i = 1:800
% fprintf(fileID,'%f %f %f\n',O2nighttext(i,1),O2nighttext(i,2)*(-1),O2nighttext(i,3));
% if rem(i,8) == 0
% fprintf(fileID,'\n');
% end
% end
% fclose(fileID);