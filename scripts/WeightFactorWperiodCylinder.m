function normWater = WeightFactorWperiodCylinder(waterfilm,x,y)
global R

%Cylindrical Y
thickness = min(waterfilm,waterfilm(2,2));

wd = zeros(7,1);
wd(3) = thickness(2,3);
wd(6) = thickness(2,1);
wd(7) = thickness(2,2);

if rem(y,2) == 0
    wd(1)  = thickness(1,1);
    wd(2) = thickness(1,2);
    wd(4) = thickness(3,2);
    wd(5) = thickness(3,1);
else
    wd(1) = thickness(1,2);
    wd(2) = thickness(1,3);
    wd(4) = thickness(3,3);
    wd(5) = thickness(3,2);    
end
center = wd(7);
newThickness = wd;
normW = sum(wd);

%Avoid shallow film
row  = find(wd<R);
if length(row) ~= 0
   for i = 1:length(row)
       center = center + wd(row(i));
       newThickness(row(i)) = 0;
   end
   newThickness(7) = center;
   normWater = newThickness/normW;   
else
    normWater = newThickness/normW;
end


end
