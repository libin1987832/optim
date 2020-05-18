clear
% n=2
% A=randn(n);
% A=A'*A;
% B=0.1*eye(n);
% C=A+B;
% xs=randn(n,1);
% xs(xs<0)=0;
% q=rand(n,1);
% qt=C*xs;
% q(xs>0)=-qt(xs>0);
% q(xs==0)=max(abs(qt))+0.1;
% save('fpis','C','xs','q','n')
% load('fpis')
n=2;
C=[1/4,1/7;1/7,1];
q=[-1/12;1/5];
xs=[1/3;0];
C*xs+q
x0=ones(n,1);
nmax=100;
for i=0:5
    [xks,ress]=splitS(C,q,1,x0,1);
    x0=xks
end


% 
% max_iter = 1;
% tol_rel  = 0.0;
% tol_abs  = 0.0;
% 
% 
% 
% for i=0:5
% [xks,ress]=splitS(C,q,1.4,x0,1);
% x0=xks
% end
% 
% [xkpsor err iter flag convergence msg] = psor(C, q, x0, 1.4, max_iter, tol_rel, tol_abs, false);
% [res2,fx2]=test_valid(C,q,xkpsor);
% [res3,fx3]=test_valid(C,q,xks);
% [ress,fxs]=test_valid(C,q,xs);