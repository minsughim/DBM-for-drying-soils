function tempMatrix = generateBCcylinderY(ObjectM)
%Periodic BC
[m n] = size(ObjectM);

tempMatrix = zeros(m+2, n+2);
%periodic for vertical direction
%Zero flux for horizontal direction

tempMatrix(2:m+1,1) = ObjectM(:,n);
tempMatrix(2:m+1, 2:n+1) = ObjectM(:,:);
tempMatrix(2:m+1, n+2) = ObjectM(:,1);
tempMatrix(1,:) = tempMatrix(2,:);
tempMatrix(m+2,:) = tempMatrix(m+1,:);

% tempMatrix(1,2:n+1) = ObjectM(m,:);
% tempMatrix(m+2,2:n+1) = ObjectM(1,:);
% tempMatrix(2:m+1,1) =  ObjectM(:,1);
% tempMatrix(2:m+1, n+2) = ObjectM(:,n);
% tempMatrix(2:m+1, 2:n+1) = ObjectM(:,:);
% 
% tempMatrix(1,1) = ObjectM(m,1);
% tempMatrix(m+2,1) = ObjectM(1,1);
% tempMatrix(m+2,n+2) = ObjectM(1,n);
% tempMatrix(1,n+2) = ObjectM(m,n);

end