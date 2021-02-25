% u = clustering labeling
% v = true labeling
function out = pairwiseindex(u,v)
u = u(:); v = v(:);
n = numel(u);
uk = max(u);
vk = max(v);
% Mismatch matrix
uh = hammdist(u');
vh = hammdist(v');
a = sum((uh==0).*(vh==0)); % In the same cluster - In the same cluster
b = sum(uh==0) - a;        % In the same cluster - In different cluster
c = sum(vh==0) - a;        % In different cluster - In the same cluster
d = sum(uh==vh) - a;       % In different cluster - In different cluster
% Contingency matrix
cm = full(sparse(u,v,1,uk,vk));
rt = sum(cm,2)';
ct = sum(cm,1);
srt = sum(rt.^2);
sct = sum(ct.^2);
nc = (n*(n^2+1)-(n+1)*srt-(n+1)*sct+2*(srt*sct)/n)/(2*(n-1));
abcd = a+b+c+d;
% Rand Index
ri = (a+d)/abcd;
% Adjusted Rand Index
if isequal(abcd,nc)
   ari = 0;
else
   ari = (a+d-nc)/(abcd-nc);
end
% Wallace coefficient
wab = a/(a+b);
wba = a/(a+c);
% Jaccard index
jrd = a/(a+b+c); 
% Fowlkes-Mallows index
fm = a/sqrt((a+b)*(a+c));
% Larsen index
sm = ct(ones(1,uk),:) + rt(ones(1,vk),:)';
lm = 2*cm./sm;
Lab = mean(max(lm,[],2));
Lba = mean(max(lm,[],1));
% Meila-Heckerman index
lm2 = lm;
K = min([uk vk]);
cval = 0;
for i = 1:K
    [rr,cc] = find(lm2==max(lm2(:)),1,'first');
    lm2(rr,:) = 0;
    lm2(:,cc) = 0;
    cval = cval + cm(rr,cc);
end
mh = cval/n;
% Mirkin coefficient
m = (2*(b+c))/(n*n);
% Salidas
% Rand index
% Adjusted rand index
% Wallace A->B
% Wallace B->A
% Jaccard
% Fowlkes-Mallows
% Larsen A->B
% Larsen B->A
% Meila-Heckerman
% Mirkin
% out = [ri,ari,wab,wba,jrd,fm,Lab,Lba,mh,m];
out = ari;
%****************************************************
function H = hammdist(A)
AdotA = dot(A,A,1);
D = bsxfun(@plus,AdotA,(AdotA)') - 2*(A'*A) > 0;
I = tril(true(size(D)),-1);
H = D(I)';