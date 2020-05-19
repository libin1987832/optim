clear
%  n=100;
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
% x0=ones(n,1);
% save('fpis','C','xs','q','n')

% load('fpis')
C=[1/4,1/7;1/7,1];
D=tril(C)
inv(D)
% x0=ones(n,1);
% xA=[];
% I0=find(x0>0);
% begin=1;
% xAI=[];
% iter=20;
% I=zeros(1,iter);
% xAI=[I];
% xAI(begin,1)=1;
% for i=2:iter
%     [xks,ress]=splitS(C,q,1.4,x0,1);
%     Iks=find(xks>0);
%     e=setdiff(I0,Iks);
%     if isempty(e)
%         xAI(begin,i)=1;
%     else
%         xAI=[xAI;zeros(1,iter)];
%         begin=begin+1;
%         xAI(begin,i)=1;
%         I0=Iks;
%     end
%     x0=xks;
%     xA=[xA x0];
% end
% xAI





% n=2;
% C=[1/4,1/7;1/7,1];
%
% xs=[1/3;1];
% q=-C*xs;
% x0=[0;0];
% nmax=100;
% xA=[];
% for i=0:5
%     [xks,ress]=splitS(C,q,1,x0,1);
%     x0=xks;
%     xA=[xA x0];
% end
% xA

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