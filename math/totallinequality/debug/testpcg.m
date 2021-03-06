% load('AF1.mat')
% A=AF1'*AF1;
% b=AF1'*b(rF1>1e-5);
% x0=sparse(zeros(130,1));
% [x0,fl0,rr0,it0,rv0] = pcg(A,b,1e-8,100);
load('test.mat')
A=full(A);
b=A'*b;
A=A'*A;

L = ichol(A);
M = L*L';
r0=b-A*x0;
z0=L\r0;
y0=L'\r0;
p0 = y0;
ry0= r0' * y0;
xa=[];
index=1
maxit = 2;
it3=maxit ;
rv3 = zeros(maxit+1,1); 
rv3(1,1)=norm(r0);
while norm(r0)>1e-8 && index < maxit;
Ap0=A*p0;  
a = ry0/(p0'*Ap0);
x1 = x0 + a * p0;
r0 = r0 - a * Ap0;
z0 = L\r0;
y0 = L'\r0;
ry1 = r0'*y0;
beta= (ry1/ry0);
p0 = y0+beta*p0;
ry0=ry1;
index = index + 1;
rv3(index,1)=norm(r0);
end
[x1p,fl1,rr1,it1,rv1] = pcg(A,b,1e-8,maxit,L,L');
[x2,fl2,rr2,it2,rv2] = pcg(A,b,1e-8,maxit,L,[]);
% L = ichol(A,struct('michol','on'));
% [x2,fl2,rr2,it2,rv2] = pcg(A,b,1e-8,100,L,L');
% 
figure;
semilogy(0:it3,rv3/norm(b),'b.');
hold on;
semilogy(0:it1,rv1/norm(b),'r.');
 semilogy(0:it2,rv2(1:it2+1)/norm(b),'k.');
legend('No Preconditioner','IC(0)','MIC(0)');
xlabel('iteration number');
ylabel('relative residual');
hold off;
