% rk = b-Ax0 rkN =£¨b-Ax0£©+ rRkN= || £¨b-Ax0£©+ ||^2 nGARk = ||
% A^T£¨b-Ax0£©+||^2
% note that computation A * x0 A' * rkN
function [rk, rkN, nRkN, nGARk] = residual(A,b,x0)
rk = b - A * x0;
rkN = rp;
rkN(rkN<0) = 0;
GARk = A' * rkN;
nGARk = GARk'*GARk;
nRkN = rkN'*rkN;
