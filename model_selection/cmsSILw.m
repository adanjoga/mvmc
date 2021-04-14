
% Clustering Model Selection (CMS) using Silhouette (Sil)
%   Clr is best clustering solution from the PF in terms of Sil
function [Clr,SILb,idx,SILvalues,Wi] = cmsSILw(CLRs, W, Dataviews, Nobj)

PFsize = size(CLRs,2);
DVnorm = cell(Nobj,1);

% This step is not needed if the data views in Dataviews are normalized
for m = 1:Nobj
    DVnorm{m} = minmaxnorm(Dataviews{m});
end

SILvalues = NaN(PFsize,1);
for i =1:PFsize
    Wi = W(i,:);
    Ci = CLRs(:,i);
    Ki = numel(unique(Ci));
    
    % Computation of the weigted dissimilarity matrix
    Dws = 0;
    for m = 1:Nobj
        Dws = Dws + (Wi(m) .* DVnorm{m});
    end
    
    % Computation of the Wil index considering Dws matrix
    SILvalues(i) = silindex_ws(Ci, Ki, Dws);
end
[SILb, idx] = max(SILvalues);
Clr = CLRs(:,idx);
Wi = W(idx,:);
end