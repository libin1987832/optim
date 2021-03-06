% clear
%  n=2;
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
% save('fpis','C','xs','q','n','x0')

% load('fpis')
% C=[4,-3;-3,8];
% q=[5;-16];
% x0=ones(n,1);
% [0,3/4;0,9/32]*[0,2]-[5,-16]
% (0.75)^2*(9/32)*2-(0.75)^2*(-49/32)-(0.75)*(-49/32)-5/4
% D=tril(C)
% dl=inv(D);
% x=dl*(D-C)*x0-dl*q
% [0,3/4;0,9/32]^3*[0;2]-[0,3/4;0,9/32]^2*[5/4;-49/32]-[0,3/4;0,9/32]*[5/4;-49/32]-[5/4;-49/32]


% xA=[];
% I0=find(x0>0);
% begin=1;
% xAI=[];
% iter=20;
% I=zeros(1,iter);
% xAI=[I];
% xAI(begin,1)=1;
% for i=2:iter
%     [xks,ress]=splitS(C,q,1,x0,1);
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
% xA

% C=[1/4,1/7;1/7,1];
% D=tril(C)
% inv(D)


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

%% predict
n=30;
cz=zeros(1,n);
for e=1:100
    A=randn(n);
    A=A'*A;
    B=0.1*eye(n);
    C=A+B;
    xs=randn(n,1)+40;
    qt=-C*xs;
    % Cx+qt=0;
    % C=[4,-3;-3,8];
    % qt=[5;-16];
    
    x0=ones(n,1);
    iter=40;
    t=predict(C,x0,qt,2);
    rx=x0;
    for i=1:10
    rx=computDLU(C,rx);
    end
    rxt=rx+t;
    t=rxt;
    xk0=x0;
    xk=x0;
    for i=1:iter
        for j=1:n
            xk(j)=xk0(j)-(C(j,:)*[xk(1:j-1);xk0(j:n)]+qt(j))/C(j,j);
        end
        xk0=xk;
    end
    I1=(xk>0);
    I2=(t>0);
    d=abs(I1-I2);
%     [ms,ns]=size(d);
    cz(sum(d)+1)=cz(sum(d)+1)+1;
end
cz