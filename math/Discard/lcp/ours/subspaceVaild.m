function [e,isposdef]=subspaceVaild(x0,M,q)
I=(x0>0);
MII=M(I,I);
qI=q(I);
x=zeros(size(x0));
e=norm(MII*x0(I)+q(I));

% tf = issymmetric(MII)
d = eig(MII);
tol=length(d)*eps(max(d));
isposdef = all(d > tol);
% issemidef = all(d > -tol)
 

