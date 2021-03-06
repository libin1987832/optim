% pina hybrid algorithm
function [xk,rk,countFM,countNW,beginNW,tf,vk]=pina(x0,A,b,maxIter)
t=clock;
tol=1e-15;
%compute hybrid uIter
[m,n]=size(A);

%FM need a qr decompose
r0=b-A*x0;
r0(r0<0)=0;
%condition for terminate
Ar=norm(A'*r0);
rn=norm(r0);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;

countFM=0;
countNW=0;
beginNW=0;


% AII=AI'*AI;
% bII=AI'*(r(I));
% sparse for svds densy for svd
% [U,S,V]=svds(AII);

if Ar<delt*rn || rn<delt
    xk=x0;
    rk=r0;
    fk=0.5*(r0'*r0);
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt
    I=find(r0>=tol);
    %提取子矩阵判断是否正定
    AI=A(I,:);
    AII=AI'*AI;
    [d,v] = eig(AII);
    isposdef = all(d > tol);
    hk=AI\r0(I);
    if isposdef
        countFM= countFM+1;
        aa=bisect(r0,A*hk);
%         aa=piecewise(A,b,dh,x)
    else
        countNW= countNW+1;
        ai=r0./(A*hk);
        aa=min(ai(ai>tol));
    end
    xk=x0+aa*hk;
    rk=b-A*xk;
    rk(rk<0)=0;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
    r0=rk;
    vk=sum(sign(rk));
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));
%disp(['%hybrid1 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);
%disp(['$',num2str(m),'\times ',num2str(n),'$&FM&(',num2str(countFM),',',num2str(countNW),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
%disp(['well1033&Daxs&',num2str(vk),'&',num2str(rn),'&',num2str(Ar),'&(',num2str(countFM),',',num2str(countNW),')&',num2str(beginNW)]);


