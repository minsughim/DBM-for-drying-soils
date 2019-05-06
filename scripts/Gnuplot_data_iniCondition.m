 % For chemical profile without cells
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

xlist = reshape(listV(:,1), [m,n])*1000;
ylist = reshape(listV(:,2), [m,n])*(-1000);
LocalGasContentstext = zeros(m*n,3);
LocalGasContentstext(:,1) = xlist(:);
LocalGasContentstext(:,2) = ylist(:);
LocalGasContentstext(:,3) = LocalGasContents(:);
LocalGasContentstext2 = sortrows(LocalGasContentstext,-2);
fileID = fopen('./../../../Fig/LocalGasContents.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',LocalGasContentstext2(i,1),LocalGasContentstext2(i,2)*(-1),LocalGasContentstext2(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);

O2daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = sitesC{1}(:,:);
O2daytext(:,1) = xlist(:);
O2daytext(:,2) = ylist(:);
O2daytext(:,3) = temp(:);
O2daytext = sortrows(O2daytext,-2);
fileID = fopen('./../../../Fig/O2profile.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


pHdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = sitesC2{5}(:,:);
pHdaytext(:,1) = xlist(:);
pHdaytext(:,2) = ylist(:);
pHdaytext(:,3) = temp(:);
pHdaytext = sortrows(pHdaytext,-2);
fileID = fopen('./../../../Fig/pHprofile.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',pHdaytext(i,1),pHdaytext(i,2)*(-1),pHdaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);

TANdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = sitesC{6}(:,:)+sitesC2{3}(:,:);
TANdaytext(:,1) = xlist(:);
TANdaytext(:,2) = ylist(:);
TANdaytext(:,3) = temp(:);
TANdaytext = sortrows(TANdaytext,-2);
fileID = fopen('./../../../Fig/TANprofile.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',TANdaytext(i,1),TANdaytext(i,2)*(-1),TANdaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


NO3daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = sitesC{5}(:,:);
NO3daytext(:,1) = xlist(:);
NO3daytext(:,2) = ylist(:);
NO3daytext(:,3) = temp(:);
NO3daytext = sortrows(NO3daytext,-2);
fileID = fopen('./../../../Fig/NO3profile.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


Islanddaytext = zeros(m*n,3);
%temp = zeros(m,n);
temp = InvasedIsland;
Islanddaytext(:,1) = xlist(:);
Islanddaytext(:,2) = ylist(:);
Islanddaytext(:,3) = temp(:);
Islanddaytext = sortrows(Islanddaytext,-2);
fileID = fopen('./../../../Fig/Islandprofile.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',Islanddaytext(i,1),Islanddaytext(i,2)*(-1),Islanddaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);



%% For chemical profile with cells at water filled capacity
%
nightT = 1368;
dayT = 1224;
HoriList = 0:Lp:Lp*(n-1);
% figure(2)
% subplot(1,2,1)
% contour(xlist,ylist,  timeConcDist{1}(:,:,dayT))
% caxis([0 35])
% colorbar;
% subplot(1,2,2)
% contour(xlist,ylist,  timeConcDist{1}(:,:,nightT))
% caxis([0 35])
% colorbar;

O2daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{1}(:,:,dayT);
O2daytext(:,1) = xlist(:);
O2daytext(:,2) = ylist(:);
O2daytext(:,3) = temp(:);
O2daytext = sortrows(O2daytext,-2);
fileID = fopen('./../../../Fig/O2profileDay.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


O2daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{1}(:,:,nightT);
O2daytext(:,1) = xlist(:);
O2daytext(:,2) = ylist(:);
O2daytext(:,3) = temp(:);
O2daytext = sortrows(O2daytext,-2);
fileID = fopen('./../../../Fig/O2profileNight.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',O2daytext(i,1),O2daytext(i,2)*(-1),O2daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);



pHdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist2{5}(:,:,dayT);
pHdaytext(:,1) = xlist(:);
pHdaytext(:,2) = ylist(:);
pHdaytext(:,3) = temp(:);
pHdaytext = sortrows(pHdaytext,-2);
fileID = fopen('./../../../Fig/pHprofileDay.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',pHdaytext(i,1),pHdaytext(i,2)*(-1),pHdaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);

pHdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist2{5}(:,:,nightT);
pHdaytext(:,1) = xlist(:);
pHdaytext(:,2) = ylist(:);
pHdaytext(:,3) = temp(:);
pHdaytext = sortrows(pHdaytext,-2);
fileID = fopen('./../../../Fig/pHprofileNight.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',pHdaytext(i,1),pHdaytext(i,2)*(-1),pHdaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


TANdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{6}(:,:,dayT)+timeConcDist2{3}(:,:,dayT);
TANdaytext(:,1) = xlist(:);
TANdaytext(:,2) = ylist(:);
TANdaytext(:,3) = temp(:);
TANdaytext = sortrows(TANdaytext,-2);
fileID = fopen('./../../../Fig/TANprofileDay.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',TANdaytext(i,1),TANdaytext(i,2)*(-1),TANdaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);

TANdaytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{6}(:,:,nightT)+timeConcDist2{3}(:,:,nightT);
TANdaytext(:,1) = xlist(:);
TANdaytext(:,2) = ylist(:);
TANdaytext(:,3) = temp(:);
TANdaytext = sortrows(TANdaytext,-2);
fileID = fopen('./../../../Fig/TANprofileNight.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',TANdaytext(i,1),TANdaytext(i,2)*(-1),TANdaytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);



NO3daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{5}(:,:,dayT);
NO3daytext(:,1) = xlist(:);
NO3daytext(:,2) = ylist(:);
NO3daytext(:,3) = temp(:);
NO3daytext = sortrows(NO3daytext,-2);
fileID = fopen('./../../../Fig/NO3profileDay.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


NO3daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{5}(:,:,nightT);
NO3daytext(:,1) = xlist(:);
NO3daytext(:,2) = ylist(:);
NO3daytext(:,3) = temp(:);
NO3daytext = sortrows(NO3daytext,-2);
fileID = fopen('./../../../Fig/NO3profileNight.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);


NO3daytext = zeros(m*n,3);
temp = zeros(m,n);
temp = timeConcDist{5}(:,:,nightT);
NO3daytext(:,1) = xlist(:);
NO3daytext(:,2) = ylist(:);
NO3daytext(:,3) = temp(:);
NO3daytext = sortrows(NO3daytext,-2);
fileID = fopen('./../../../Fig/NO3profileNight.txt','w');
for i = 1:m*n
fprintf(fileID,'%f %f %f\n',NO3daytext(i,1),NO3daytext(i,2)*(-1),NO3daytext(i,3));
if rem(i,20) == 0
fprintf(fileID,'\n');
end
end
fclose(fileID);

