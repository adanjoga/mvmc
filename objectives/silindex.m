function f = silindex(Clrs,K,DXX,lowerb)
% Validacion del agrupamiento
N = numel(Clrs);
clusts = unique(Clrs);
Nk = accumarray(Clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<3)
    f = lowerb;
    return;    
end

% Evaluacion del IVG
DX = DXX.^2;
S = zeros(1,K);
for i = 1:K
    XC = DX(Clrs==clusts(i),Clrs==clusts(i));
    ai = sum(XC,2)/max(Nk(i)-1, 1); % modificacion 25-10-16
    %ai = mean(XC,2);
    dd = zeros(Nk(i),K-1);
    k = 0;
    for j = setdiff(1:K,i)
        k = k+1;
        XK = DX(Clrs==clusts(i),Clrs==clusts(j));
        dd(:,k) = mean(XK,2);
    end
    bi = min(dd,[],2);
    S(i) = sum((bi-ai)./max([ai bi],[],2),1);
    
    % Restriccion:  The score is 0 for clusters with size = 1
    % Modificacion: 26 de Noviembre de 2018
    if Nk(i) == 1 
        S(i) = 0;
    end
end
f = (1/N)*sum(S);


