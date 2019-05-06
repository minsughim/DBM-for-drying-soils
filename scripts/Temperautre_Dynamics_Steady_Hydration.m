DeltaT = 20000;
dampD =sqrt(mean(alphaMatrix(:))*3600*24/3.14);

testMinTt = 0.5*Lp^2/(4*max(alphaMatrix(:)));

diffMatrixInd = alphaMatrix;
DiagDiff = diag(diffMatrixInd(:));
nnDiff = listE*DiagDiff;
for i = 1:(m*n)
    nnEffDiff(:,i) = 2*nnDiff(:,i).*diffMatrixInd(i)./(nnDiff(:,i)+diffMatrixInd(i));
end
tempN = DeltaT/testMinTt;
tempSg = zeros(m,n);

%sitesT = 15*ones(m,n); %global soil temperature
for i = 1:n
    sitesT(:,i) =  averageT + amplitudeT.*exp(-depthList./dampD).*sin((2*pi*(timeLine(1)-6)/24)-(depthList./dampD));
end
tempST = sitesT;
HeteroDiffM = nnEffDiff*testMinTt/Lp/Lp; %nearest neighbor
HeteroDiffM = HeteroDiffM - diag(sum(HeteroDiffM)); %including self
Afuture = sparse((eye(m*n)-HeteroDiffM));
Apresent = (eye(m*n)+HeteroDiffM);
for i = 1:tempN
    tempST(1,:) = avgT;
    tempST(end,:) = avgT;
    presentChange = Apresent*tempST(:);
    [tempC1,flag] = bicgstab(Afuture,presentChange, [],1);
    tempST = reshape(tempC1,m,n);
end
sitesTSteady = reshape(tempST,m,n);
%hold on
%mesh(sitesSteady)

DeltaT = dt*plottt;
tempN = DeltaT/testMinTt;
sitesTevolove = zeros(m,n,dataPoints);
for timeT = 1:dataPoints
    %timeT
    for i = 1:tempN
        tTemp = timeLine(timeT);
        tempST(1,:) = averageT + amplitudeT.*sin(2*pi*(tTemp-6)/24);
        tempST(end,:) = averageT + amplitudeT.*exp(-depthList(end)./dampD).*sin((2*pi*(tTemp-6)/24)-(depthList(end)./dampD));
        presentChange = Apresent*tempST(:);
        [tempC1,flag] = bicgstab(Afuture,presentChange, [],1);
        tempST = reshape(tempC1,m,n);
    end
    sitesTevolove(:,:,timeT) = reshape(tempST,m,n);
end
% 
% meanTempdist = zeros(m,dataPoints); 
% for timeT = 1:dataPoints
% index = index + 1;
% temp = sitesTevolove(:,:,timeT);
% meanTempdist(:,index) = mean(temp,2);
% end
% 
