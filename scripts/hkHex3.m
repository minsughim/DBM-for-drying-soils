function [result] = hkHex3(input)

% Apply the Hoshen-Kopelman Algorithm for Cluster Identification
matrix = round(input);

[m,n] = size(matrix);

% labels = zeros(m*n/2,1);  % allocate space for the labels

n_labels = 0;             % number of labels used so far
labels=[];
removedL = [];
%  /* scan the matrix */

for i=1:m,
    for j=1:n
        if matrix(i,j) == 1
            
            % check the neighbors of this cell, up and to the left
            up = 0;
            left = 0;
            ne = 0;            
            if (i > 1),
                up = matrix(i-1,j);
            end
            if (j > 1)
                left = matrix(i,j-1);
            end
            if rem(i,2) ~= 0
                if (i > 1)&&(j < n)
                    ne = matrix(i-1,j+1);
                end
            else
                if (j > 1)&&(i > 1)
                    ne = matrix(i-1,j-1);
                end
            end                       
            temp = ((up>0)  + (left>0) + (ne>0));            
            if temp == 0
                n_labels = n_labels + 1;           % make the new label number
                labels( n_labels ) = n_labels;     % no known alias yet (label[x]=x)
                matrix(i,j) = n_labels;            % apply this label to this position
            else
                matrix(i,j) = max([up,left,ne]);      % whichever is nonzero is already labelled; use that label                
            end % switch(up+left)
        end % if
    end % for
end % for

%  Renumber the labels so that they're continuous, and eliminate cluster aliases

result1 = matrix;

for i=1:m,
    for j=1:n
        
        selfL = matrix(i,j);
        
        if selfL > 0
            
            % check the neighbors of this cell, up and to the left
            down = 0;
            right = 0;
            sw = 0;
            up = 0;
            left = 0;
            ne = 0;
            
            if (i > 1),
                up = matrix(i-1,j);
            end
            if (j > 1)
                left = matrix(i,j-1);
            end
            
            if rem(i,2) ~= 0
                if (i > 1)&&(j < n)
                    ne = matrix(i-1,j+1);
                end
                if (i < m)&&(j < n)
                    sw = matrix(i+1,j+1);
                end
                
            else
                if (j > 1)&&(i > 1)
                    ne = matrix(i-1,j-1);
                end
                if (j > 1)&&(i < m)
                    sw = matrix(i+1,j-1);
                end
                
            end
            
            if (i < m),
                down = matrix(i+1,j);
            end
            if (j < n)
                right = matrix(i,j+1);
            end
            
            temp = ((down>0)  + (right>0) + (sw>0) + (up>0)  + (left>0) + (ne>0));
            
            if temp ~= 0
                
                tempLabels =[up,down,left,right,sw,ne,selfL];
                [c,Ilist] = sort(tempLabels,2);
                list = find(c>0);
                label = Ilist(list(1));
                tempMin = tempLabels(label);
                for ik = 1:7
                    if (tempLabels(ik)~=0)&&(tempLabels(ik)~=tempMin)
                        tempMax = tempLabels(ik);
                        matrix(find(matrix == tempMax)) = tempMin;
                    end
                end
                
            end % switch(up+left)
        end % if
    end % for
end % for

result2 = matrix;

maxLable = max(max(matrix));
a = 1;
for iL = 1:maxLable
    [tempC,tempI] = find(matrix == iL);
    if length(tempC) == 0
        removedL(a) = iL;
        a = a+1;
    end
end

correctionL = length(removedL);

if correctionL ~= 0
    for iL1 = 1:correctionL
        matrix(find(matrix == (maxLable+1-iL1))) = removedL(correctionL+1-iL1);
    end
end

result = matrix;


% %disp('relabeling'),
% j = 1;
% for i=1:n_labels
%     if (labels(i) == i)
%         labels(i) = j;
%         j = j+1;
%     else
%         labels(i) = labels(labels(i));
%     end
% end
%
% %  apply the relabeling to the matrix
%
% for i=1:m
%     for j=1:n
%         if not (matrix(i,j) == 0)
%             matrix(i,j) = labels(matrix(i,j)+0);
%         end
%     end
% end

%result = matrix;

