function D = eucdist(A, B)
% D = sqrt(bsxfun(@plus,dot(B,B,1),dot(A,A,1)')- 2*mtimesx(A,'T',B));
D = sqrt(bsxfun(@plus,dot(B,B,1),dot(A,A,1)')- 2*(A'*B));