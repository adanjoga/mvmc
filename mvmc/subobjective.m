function obj = subobjective(weight, ind, idealpoint, method)
%SUBOBJECTIVE function evaluate a point's objective with a given method of
%decomposition. 

%   Two method are implemented by far is Weighted-Sum and Tchebesheff.
%   weight: is the decomposition weight.(column wise vector).
%   ind: is the individual point(column wise vector).
%   idealpoint: the idealpoint for Tchebesheff decomposition.
%   method: is the decomposition method, the default is 'te' when is
%   omitted.
%   
%   weight and ind can also be matrix. in which have two scenairos:
%   When weight is a matrix, then it's treated as a column wise set of
%   weights. in that case, if ind is a size 1 column vector, then the
%   subobjective is computed with every weight and the ind; if ind is also
%   a matrix of the same size as weight, then the subobjective is computed
%   in a column-to-column, with each column of weight computed against the
%   corresponding column of ind. 
%   A row vector of subobjective is return in both case.

    if (nargin==2)
        obj = ws(weight, ind);
    elseif (nargin==3)
        obj = te(weight, ind, idealpoint);
    else
        if strcmp(method, 'wsu')
            obj=ws(weight, ind);
        elseif strcmp(method, 'tch')
            obj=te(weight, ind, idealpoint);
        elseif strcmp(method, 'pbi_mn')
            obj=pbi_mn(weight, ind, idealpoint);
        elseif strcmp(method, 'pbi_mx')
            obj=pbi_mx(weight, ind, idealpoint);
        else
            obj= te(weight, ind, idealpoint);
        end
    end
end

function obj = ws(weight, ind)
    if size(ind, 2) == 1 
       obj = (weight'*ind)';
    else
       obj = sum(weight.*ind);
    end
end

function obj = te(weight, ind, idealpoint)
    s = size(weight, 2);
    indsize = size(ind,2);
    
    weight((weight == 0))=0.00001;
    
    if indsize==s 
        part2 = abs(ind-idealpoint(:,ones(1, indsize)));
        obj = max(weight.*part2);
    elseif indsize ==1
        part2 = abs(ind-idealpoint);
        obj = max(weight.*part2(:,ones(1, s)));   
    else
        error('individual size must be same as weight size, or equals 1');
    end
end

%% Penalty-based boundary intersection (PBI) function for MINIMIZATION problems
function obj = pbi_mn(weight, ind, idealpoint)
    s = size(weight, 2);
    indsize = size(ind,2);
    
    weight((weight == 0))=0.00001;
    sita = 0.9; % Parameter value suggested by the author pp. 724 (left)
    
    if indsize==s
        dd1 = (ind-idealpoint(:,ones(1,indsize))) .* weight;
        d1 = sqrt(sum(dd1.*dd1,1)) ./ (sqrt(sum(weight.*weight,1)));
        %d1=abs(sum((ind-idealpoint(:,ones(1,indsize))).*weight,1)) ./ (sqrt(sum(weight.*weight,1))); 
        dd2=ind-(idealpoint(:,ones(1,s)) + d1(ones(size(weight,1),1),:).*weight);
        d2=sqrt(sum(dd2.*dd2,1));
        obj=d1 + sita * d2;
    elseif indsize==1
        d1t = ind-idealpoint;
        dd1 = d1t(:,ones(1,s)).* weight;
        d1 = sqrt(sum(dd1.*dd1,1)) ./ (sqrt(sum(weight.*weight,1)));
        % d1=abs(sum(d1t(:,ones(1,s)).*weight,1)) ./ sqrt(sum(weight.*weight,1));
        dd2=ind(:,ones(1,s)) - (idealpoint(:,ones(1,s)) + d1(ones(size(weight,1),1),:).*weight);
        d2=sqrt(sum(dd2.*dd2,1)); % Euclidean Norm
        obj = d1 + sita * d2;
    else
        error('individual size must be same as weight size, or equals 1');
    end
end


%% Penalty-based boundary intersection (PBI) function for MAXIMIZATION problems
function obj = pbi_mx(weight, ind, idealpoint)
    s = size(weight, 2);
    indsize = size(ind,2);
    
    weight((weight == 0))=0.00001;
    sita = 5; % Parameter value suggested by the author pp. 724 (left)
    
    if indsize==s
        d1=abs(sum((idealpoint(:,ones(1,indsize))-ind).*weight,1)) ./ (sqrt(sum(weight.*weight,1))); 
        dd=ind-(idealpoint(:,ones(1,s)) - d1(ones(size(weight,1),1),:).*weight);
        d2=sqrt(sum(dd.*dd,1));
        obj=d1 + sita * d2;
    elseif indsize==1
        part2 = abs(idealpoint-ind);
        d1=abs(sum(part2(:,ones(1,s)).*weight,1)) ./ sqrt(sum(weight.*weight,1));
        dd=ind(:,ones(1,s)) - (idealpoint(:,ones(1,s)) - d1(ones(size(weight,1),1),:).*weight);
        d2=sqrt(sum(dd.*dd,1)); % Euclidean Norm
        obj = d1 + sita * d2;
    else
        error('individual size must be same as weight size, or equals 1');
    end
end