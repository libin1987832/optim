function [xk,rk,countFM,countNW,beginNW,tf,vk]=residualR(x0,A,b,maxIter)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
%FM need a qr decompose
r0=b-A*x0;
rkp=r0;
r0(r0<0)=0;

[Q,R]=qr(A);
Qn=Q(:,1:n);
%condition for terminate
Ar=norm(A'*r0);
rn=norm(r0);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
uIndex=0;

countFM=0;
countNW=0;
beginNW=0;

if Ar<delt*rn || rn<delt
    xk=x0;
    rk=r;
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end

%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    Qr=Qn'*r0;
    rkp=rkp-Qn*Qr;
    rk=rkp;
    rk(rk<0)=0;
    r0=rk;
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    countFM=countFM+1;
    if maxIter < countFM
        break;
    end
end
xk=A\(b-rkp);
tf=etime(clock,t);
vk=sum(sign(rk));