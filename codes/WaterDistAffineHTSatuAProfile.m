function [waterFilm,percolProb,waterSatu,totalPoreVolM,sepecifInterA,LengthList] = WaterDistAffineHTSatuAProfile(systemAffine,pot,H,R1,R2,R)
%Calcualtion with the domain with constant D :: specific interfacial area
%included
%global R1 R2 R

[m,n,l] =size(systemAffine);


% %distribution of percolation threshold (C.Du, C. Satik, and Y.C. Yotsos
% %1996 : Percolation in a Fractional Browian Motion lattice)
% mu1 = (H==1)*0.386 + (H==0.8)*0.39 + (H==0.65)*0.405 + (H==0.5)*0.42 + (H==0.35)*0.44 + (H==0.2)*0.46 + (H==0)*0.5;
% sigma1 = (H==1)*0.13 + (H==0.8)*0.12 + (H==0.65)*0.115+(H==0.5)*0.11+(H==0.35)*0.1+(H==0.2)*0.08+(H==0)*0;
% %pd = makedist('Normal','mu', mu1,'sigma',sigma1);

constFilm = -1*1.9*10^(-19)/(9.81*6*pi*998);
vdwFilmcontri = (constFilm/pot)^(1/3);
r = 7.41*10^(-6)/(-pot);
crit = 2*(r+vdwFilmcontri);
a = (3*R)-(6*vdwFilmcontri);
b = -3*(4-pi)*r^2;
c = -1*crit^3 +3*(4-pi)*r^2*crit + 6*vdwFilmcontri*crit^2;
P = [a b c];
rR = roots(P);
r1 = rR(1);
if r1 > R1
    r1 = R1;
elseif r1 < 3*R
    r1 = 3*R;
end
waterFilm = zeros(m,n);
percolProb = zeros(m,n);
waterSatu = zeros(m,n);
totalPoreVolM = zeros(m,n);
sepecifInterA = zeros(m,n);
x=sym('x');
D1 = systemAffine(1,1,2);
tempA = x^(2-D1);
FilmA = int(tempA, x, R2, R1);
FilmAbsorbed = vdwFilmcontri*double(FilmA);

portionA = int(tempA, x, 3*R, r1);
ListNa = (x)^(-D1);

indiSaturation = (crit/x)^3 + 3*(4-pi)*((r/x)^2)*(1-crit/x)+6*(vdwFilmcontri/x)*(1-(crit/x)^2);
totalV = ((x^3)/3)*ListNa;
unsatuV = ((x^3)/3)*indiSaturation*ListNa;
satuArea = int(totalV, x, R2, crit);
unsatuArea = int(unsatuV, x, crit, R1);
totalPoreVol = int(totalV, x, R2, R1);
TotalVolume = double(satuArea) + double(unsatuArea);

%Expected tortuousity
totalL = x*ListNa;
waterL = x*indiSaturation*ListNa;
totaltravelL = int(totalL, x, R2, crit);
totaltortuousL = int(waterL, x, crit, R1);
LengthList(1) = double(totaltravelL);
LengthList(2) = double(totaltortuousL);

InterfacialA = 2*(x-crit)^2 + 2*pi*r*(x-crit)+ 2*pi*r*((1-sqrt(5))*r+crit);
totalIntfA = (x^2)*ListNa;
unsatuIntfA = InterfacialA*ListNa;
satuIntfArea = int(totalIntfA, x, R2, crit);
unsatuIntfAreaPore = int(unsatuIntfA, x, crit, R1);
unsatuIntfAreaSolid = int(totalIntfA, x, crit, R1);
IntfAreaSec(1) = double(satuIntfArea);
IntfAreaSec(2) = double(unsatuIntfAreaPore);
IntfAreaSec(3) = double(unsatuIntfAreaSolid);

if double(satuArea) > double(totalPoreVol)
    totalPoreV1 = double(totalPoreVol);
    IntfAreaSec(1) = double(FilmA);
    IntfAreaSec(3) = 0;
    IntfAreaSec(2) = 0;
    
else
    totalPoreV1 = TotalVolume;
end

FilmThickTemp = totalPoreV1/double(FilmA);

if r1 < 3*R
    partialSatuTemp = 0;
else
    partialSatuTemp = double(portionA)./FilmA;
end


for i = 1:m
    for j = 1:n
                
        beta = systemAffine(i,j,1);
        FilmThick = beta*FilmThickTemp + (1-beta)*vdwFilmcontri;
        waterFilm(i,j) = double(FilmThick);
        partialSatu = beta*partialSatuTemp;
        temp = systemAffine(i,j,3);
        if temp > partialSatu
            percolProb(i,j) = partialSatu;
        else
            percolProb(i,j) = 1;
        end
        saturation = (beta*totalPoreV1+(1-beta)*FilmAbsorbed)/(beta*double(totalPoreVol)+ (1-beta)*FilmAbsorbed);
        waterSatu(i,j) = double(saturation);
        totalPoreVolM(i,j) = ((beta*double(totalPoreVol)+ (1-beta)*FilmAbsorbed)/double(FilmA));
        
        interA = (IntfAreaSec(1) + beta*IntfAreaSec(2)+ (1-beta)*IntfAreaSec(3));
        
        if double(saturation) < 1
            sepecifInterA(i,j) = beta*interA/(totalPoreVolM(i,j)*3*double(FilmA));
        else
            sepecifInterA(i,j) = 0;
        end
        
    end
end


end
