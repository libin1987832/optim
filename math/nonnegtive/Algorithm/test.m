% A=diag([10^15,1]);
% % tf = issymmetric(A);
% 
% d = eig(A);
% tol=length(d)*eps(max(d))
% isposdef = all(d > tol)
% issemidef = all(d > -tol)
% 
% 
n=5;
A=rand(n+2,n);
[q,r]=qr(A);
Q1=q(:,1:3);
N=diag([1,1,0,0,1,1,1])
diag(ones(7,1))-Q1*Q1'*N