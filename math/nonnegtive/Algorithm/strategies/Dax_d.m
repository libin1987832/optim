% Dax hybrid algorithm
function [xk,rk,countFM,countNW,beginNW,tf,vk]=Dax(x0,A,b,maxIter)
t=clock;

%compute hybrid uIter
[m,n]=size(A);
uIter=floor(max(33,(m+n)/4));
%FM need a qr decompose
[Q,R]=qr(A);
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
    fk=0.5*(r'*r);
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    if uIndex<uIter
        countFM=countFM+1;
        r0=rkp;
        %FM algorithm
        [xk,rkp]=FixedM(x0,Q,R,A,b,r0);
        rk=rkp;
        AA=find(rkp>0);
        rk(rk<0)=0;
    else
        countNW=countNW+1;
        if countNW ==1
            beginNW=countFM;
        end
        uIndex=0;
        [xk,~]=krylov(A,b,xk,rkp);
        rk=(b-A*xk);
        rk(rk<0)=0;
    end
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));
%disp(['%hybrid1 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);
%disp(['$',num2str(m),'\times ',num2str(n),'$&FM&(',num2str(countFM),',',num2str(countNW),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
%disp(['well1033&Daxs&',num2str(vk),'&',num2str(rn),'&',num2str(Ar),'&(',num2str(countFM),',',num2str(countNW),')&',num2str(beginNW)]);


