%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain number of actors in each cluster
% input : Hoshen-Kopelman matrix
% output : n_c for every clusters
%          positive integers mean clusters of +1
%          negative integers mean clusters of -1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function results = StatisticClusters(hkTown);
% 
% temp = min(hkTown);
% numberOfMinusClusters = abs(min(temp));

temp = max(hkTown);
numberOfPlusClusters = max(temp);

results = zeros(numberOfPlusClusters,2);

for i = 1:numberOfPlusClusters

    temp = find(hkTown == i);
    results(i,2) = length(temp);
    results(i,1) = i;
    
end

% for i = 1:numberOfMinusClusters
%  
%      temp = find(hkTown == (-1)*i);
%      results(i+numberOfPlusClusters,2) = length(temp);
%      results(i+numberOfPlusClusters,1) = (-1)*i;
%      
%  end

end




