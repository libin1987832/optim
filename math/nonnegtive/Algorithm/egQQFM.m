A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
% QR decomposition in the exact arithmetic


% x=[-3;-1];
x=[-3;1]/4;
r0=b-A*x;
r0A=r0;
r0B=r0;
r0A(r0A>0)=1;
r0A(r0A<0)=0;
% violate the constrain
NK=diag(r0A)
% r_{k+1}=(I-QQ^{T}Nk)r_{k}
IQ=diag([1,1,1,1])-Q*Q'*NK;

% matlab in inexact arithmetic
[q,r]=qr(A);
for i=1:10
[xk,r0,rk,fkFM,fm,fr]=FM(x,q,r,A,b);
x=xk;
end
% % r_{k}
% r0IQ=b-A*x;
% % r_{k+1}
% rkIQ=IQ*r0IQ;
% rkIQS=rkIQ;
% rkIQS(rkIQ>0)=1;
% rkIQS(rkIQ<0)=0;
% % N_{k+1}
% NK2=diag(rkIQS)
% r_{k+1}=b-Ax_{k+1}
rkIQ2=b-A*xk;