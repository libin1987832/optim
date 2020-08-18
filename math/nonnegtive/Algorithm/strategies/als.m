function [xk,rk,countFM,countNW,beginNW,tf,vk,Arr,rkrr]=als(x0,A,b,maxIter)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
%FM need a qr decompose
r0=b-A*x0;
rkp=r0;
r0(r0<0)=0;
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
zk0=-rkp;
zk0(zk0<0)=0;
vk0=sum(sign(r0));
Arr=[];
rkrr=[];
while Ar>delt*rn && rn>delt
    xk=A\(b+zk0);
    rkp=b-A*xk;
    rkrr=[rkrr rkp];
    
    zk0=-rkp;
    rk=rkp;
    rk(rk<0)=0;
    zk0(zk0<0)=0;
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    countFM=countFM+1;
    vk=sum(sign(rk));
    if vk~=vk0
        Arr=[Arr [countFM;vk]];
        vk0=vk;
    end
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));