% m=2000;
% n=0.1*m;
% A=rand(m,n)*2-ones(m,n);
% b=rand(m,1)*2-ones(m,1);
% x0=zeros(n,1);
%
% [xk1,fk1,xkArr1]=hybrid1(x0,A,b);
% [xk2,fk2,xkArr2]=hybrid2(x0,A,b);
%[xk3,fk3,xkArr3]=hybrid3(x0,A,b);


% A=[-1,-1;1,1;-1,0;1,0;0,1;0,-1];
% b=[1,1,0.5,-0.5,0.5,-0.5]'
% x0=zeros(2,1);

% x0=-10;
% A=[1;-1];
% b=[100;2];
% [Q,R]=qr(A);
% [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b)
% m=4;
% n=2;
% A=[1,1;-1,-1;-1,0;-6,-3];
% b=[1;1;0.5;2];
% x0=[-57/100;47/100];
%  x0=[-3/4;-1/4];
% x0=[-3/4;-100];
% [Q,R]=qr(A);
% [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
% r=b-A*xk;
% Nk=r;
% Nk(Nk>0)=1;
% Nk(Nk<0)=0;
% tt.mat have feasible solution but the change very ...
% ddf.mat test result is good
% addpath('FM')
% m=10;
% ratio=0.3;
% n=ceil(ratio*m);
% A=2*rand(m,n)-1;
% b=2*rand(m,1)-1;
addpath('FM');
addpath('util');
[A,rows,cols,entries,rep,field,symm]=mmread('./util/well1033.mtx');
n=cols;
m=rows;
x0=zeros(n,1);
b=ones(rows,1);
x0=zeros(cols,1)+12441233;
save('well1033.mat','A','b','m','n')
% load('tt1.mat')
% load('ddf.mat');
[Q,R]=qr(A);
% x0=zeros(n,1);
r=b-A*x0;
Nk=r;
Nk(Nk>0)=1;
Nk(Nk<0)=0;
out=[];
fkO=[];
active=[];
QB=Q(:,1:n);
% x10=x0;
% n=10;
% B=diag(ones(m,1))-QB*QB'*diag(Nk);
%         B33=B^n;
%        B3=diag(ones(m,1))-n*QB*QB'*diag(Nk);
%        Bd=B33-B3;
%         norm(Bd*r)
%        return;
for i = 1:500
    %[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
    [xk,r0,rk,z,fk,fr,fu,fz,fd]=FMTD(x0,Q,R,A,b);
    
    if mod(i,10)==0
        % fu/fk
        rk=b-A*xk;
        % groundtrue
        r10=b-A*x10;
        
        Nk=r10;
        Nk(Nk>0)=1;
        Nk(Nk<0)=0;
        groundtrue1=sum(sign(r10.*rk));
        
        % if exposed face there is
        B=diag(ones(m,1))-QB*QB'*diag(Nk);
        B33=B^10;
        r33=B33*r10;
        groundtrue2=sum(sign(r33.*r10));
        
        % the approxim
        B3=diag(ones(m,1))-10*QB*QB'*diag(Nk);
        r3=B3*r10;
        groundtrue3=sum(sign(r3.*r10));
        disp(['i:',num2str(i),'   g1:',num2str(groundtrue1),'   r1:',num2str(groundtrue2),'  a:',num2str(groundtrue3)])
        x10=xk;
    end
    % [lambda,s]=eig(B);
    % s=diag(s);
    % less1=logical(s>1+eps*100);
    % less2=logical(s<-eps*100);
    % if sum(less1)>0
    %     error("lambda grt 1");
    % end
    % if sum(less2)>0
    %     error("lambda grt 0 ");
    % end
    
    rk=b-A*x0;
    Nk=rk;
    Nk(Nk>0)=1;
    Nk(Nk<0)=0;
%     groundtrue1=sum(sign(r10.*rk));
    r=b-A*xk;
    Nk2=r;
    Nk2(Nk2>0)=1;
    Nk2(Nk2<0)=0;
    fkO=[fkO,[i;fk;sum(Nk2)]];
    
    nkk=logical(Nk2==Nk);
    if sum(nkk) < m
        
        out=[out,[i;fk]];
        %         break;
    end
    x0=xk;
    %     activef=find(Nk2>0)';
    active=[active;Nk2'];
end
active;
out;
fkO;
base=out(1,end);
save('baseActivewell.mat','active','base');
for i=base:99
    ss=sum(active(base,:)-active(i,:));
    if ss~=0
        ss
    end
end
% r
% diag(ones(4,1))-Q(:,[1,2])*Q(:,[1,2])'*diag([1,1,0,0])
% for i = 1:20
%     [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
% x0=xk;
% rk
% end
% x0 = [-11/20;11/50];
% A = [1,1;-1,-1;-1,0;-8,-3];
% b = [1;1;1/2;3];
% [xk,fk,y]=ssqr2(x0,A,b)
% b-A*xk
% [xk,fk,y]=ssqr2(xk,A,b)
% b-A*xk
% [xk,fk,y]=ssqr2(xk,A,b)
% b-A*xk
% [xk,fk,y]=ssqr2(xk,A,b)
% b-A*xk



% A = [-1,-1;1,1;-1,0];
% b = [1,1,0.5]';
% x0 = [-10,-10]';
%
% ATA = A'*A;
% [v,d]=eig(ATA);
% lambal = max(max(d));
% for i = 1:30
%     [xk,r0,rk,fk,fm,fr] = BFM(x0,lambal,A,b);
%     x0 = xk
% end
% [xk1,fk1,xkArr1]=hybrid1(x0,A,b);
% [xk2,fk2,xkArr2]=hybrid2(x0,A,b);
% [xk3,fk3,xkArr3]=hybrid3(x0,A,b);