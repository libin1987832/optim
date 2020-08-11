function [xk,rk,countFM,countNW,beginNW,tf,vk]=Lei(x0,A,b,maxIter)
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
k=3;
countFM=0;
countNW=0;
beginNW=0;

if Ar<delt*rn || rn<delt
    xk=x0;
    rk=r0;
    fk=0.5*(r0'*r0);
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    [uk,fk]=krylov(A,rk,k);
    xk=x0+uk;
    rk=b-A*xk;
    rk(rk<0)=0;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
    countFM=countFM+1;
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));