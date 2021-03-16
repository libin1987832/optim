A=sparse([2,2 4;2,3, 5;4 5 10]);
b=sparse([2;4;5]);
x0=[0;0;0];
% rng default
% % A = sprand(400,400,.5);
% % A = A'*A;
% % b = sum(A,2);
A = delsq(numgrid('S',30));
b = ones(size(A,1),1);
x00 = zeros(size(A,1),1);
tol=1e-5;
maxit=10;
% [x,fl0,rr0,it0,rv0] = pcg(A,b,tol,maxit);
% L = ichol(A);
% [x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,maxit,L,L');
% H=inv(L);
% [x3,f3,rr3,it3,rv3] = pcg(H*A*H',H*b,tol,maxit);
% x3=H'*x3;
% rv3
L = ichol(A,struct('michol','on'));
[x2,fl2,rr2,it2,rv2] = pcg(A,b,1e-2,3,L,L',x00);
% semilogy(0:length(rv0)-1,rv0/norm(b),'-o')
% hold on
% semilogy(0:length(rv1)-1,rv1/norm(b),'-o')
% semilogy(0:length(rv2)-1,rv2/norm(b),'-o')
% % yline(tol,'r--');
% legend('No Preconditioner','Default ICHOL','Modified ICHOL','Tolerance','Location','East')
% xlabel('Iteration number')
% ylabel('Relative residual')\
x0=x00;
M=L*L';
r0=A*x0-b;
y0=L\r0;
y0=L'\y0;
p0=-y0;
ry0=r0'*y0;
res=zeros(10,1);
for i=1:10
res(i)=norm(r0);
Ap=A*p0;
a0=ry0/(p0'*Ap);
x0=x0+a0*p0;
r0=r0+a0*Ap;
y0=L\r0;
y0=L'\y0;
ry1=r0'*y0;
b0=ry1/ry0;
p0=-r0+b0*p0;
ry0=ry1;
end
