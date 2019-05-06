function [realIsland,town2,results2,numberOfCluster] = IslandStatHex(waterfilm,threshold)
[m,n] = size(waterfilm);
tempList = find(waterfilm>threshold);
town2 = zeros(m,n);
town2(tempList) = 1;
for i =1:length(tempList)
    temp = tempList(i);
    town2(temp) = 1;
end
islands2 = hkHex3(town2);
realIsland = islands2.*town2;
results2 = StatisticClusters(islands2);
listCluster = find(results2(:,2)>0);
numberOfCluster = length(listCluster);
end