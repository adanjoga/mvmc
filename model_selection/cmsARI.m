% Clustering Model Selection (CMS) using ARI
%   Clr is the best clustering solution from the PF in terms of ARI
function [Clr,ARIb,idx,ARIvalues] = cmsARI(CLRs, TLabels)

PFsize = size(CLRs,2);
ARIvalues = zeros(PFsize,1);
for i = 1:PFsize
    Yb = CLRs(:,i);
    ARIvalues(i) = pairwiseindex(TLabels,Yb);
end

[ARIb,idx] = max(ARIvalues,[],1);
Clr = CLRs(:,idx);
end