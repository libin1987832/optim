%function primal_dual(x,y)
clc;clear;

dim=4;
maxit=20;

A=[8 -1 0 0; -1 6 -3 0; 0 -3 4 -2; 0 0 -2 5];
x0=[1 0 1 0]';

c=0.1;
lowbound=0;
B = diag(ones(dim,1));

[x,crit, iters, nitBB, error] = SPL(A, B, x0, 100,  1e-12, 1e-6, 1e-9, 0, 0);
lambda=x'*A*x/(x'*B*x);
x
M = A-lambda*B;
y = M*x0
%x=Xproject(x,lowbound);

% for k=1:maxit
%     dx=zeros(dim,1);
%     x_old=x;
%     lambda=x'*A*x/(x'*x);
%     M = A-lambda*B
%     y = M*x;
%     [act,nact]=activeset(x,y,c)
%     nactMatrix = M(nact==1,nact==1)
%     %% solve
%     x(act==1)=lowbound;
%     dx_nact=-inv(nactMatrix)*y(nact==1)
%     dx(find(nact==1))=dx_nact;
%     x=x+dx;
% %    x=Xproject(x,lowbound);
%     k=k+1;
%     normsol = norm(x-x_old,2);
%     if(normsol<1e-6)
%        break; 
%     end
%     
% end
% 
% k
% lambda
% x
% y
% 
% x'*y

