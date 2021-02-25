function f = wgsindex(Clrs,K,DXX,wvalue)
% Validacion del agrupamiento
N = numel(Clrs);
clusts = unique(Clrs);
Nk = accumarray(Clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<3)
    f = wvalue;
    return;
end

sumD = zeros(K,1);
for i = 1:K
    XC = DXX(Clrs==clusts(i),Clrs==clusts(i));
    ai = sum(XC,2);
    sumD(i) = min(ai,[],1); 
    %sumD(i) = median(ai,1);
    %sumD(i) = mean(ai,1); 
end
f = sum(sumD);

