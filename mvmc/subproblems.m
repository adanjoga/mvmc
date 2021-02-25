function subp = subproblems(popsize, niche, w)
% init_weights function initialize a pupulation of subproblems structure
% with the generated decomposition weight and the neighbourhood
% relationship.
    subp = [];
    for i = 1:popsize
        p = struct('weight',[],'neighborhood',[],'individual', []);
        p.weight = w(i,:)';
        subp = [subp p];
    end
    
    %Set up the neighbourhood.
    leng=length(subp);
    distanceMatrix=zeros(leng, leng);
    for i=1:leng
        for j=i+1:leng
            A=subp(i).weight;B=subp(j).weight;
            distanceMatrix(i,j)=(A-B)'*(A-B);
            distanceMatrix(j,i)=distanceMatrix(i,j);
        end
        [~,sindex]=sort(distanceMatrix(i,:));
        subp(i).neighborhood=sindex(1:niche)';
    end
   
end