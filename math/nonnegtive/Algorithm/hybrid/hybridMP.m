function [xk,fk,xkArr,countFM,countNW,Q]=hybrid1(x0,A,b,maxIter)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
uIter=floor(max(33,(m+n)/4));
%FM need a qr decompose
[Q,R]=qr(A);
r=b-A*x0;
r(r<0)=0;
%condition for terminate
Ar=norm(A'*r);
rn=norm(r);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
uIndex=0;
xkArr=[];

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
    f1x=(x0-l)./a;
    fxc=find(f1x>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx
       sqrt()
       disp(['conject gradient:',num2str(r'*r)]);
    else
        countFM=countFM+1;
        %FM algorithm
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
        disp(['go to implement space:',num2str(r'*r)]);
    end
    k=k+1;
    bxn=bx'*bx;
    x0=xk;
end
    
    
    
    if uIndex<uIter
        countFM=countFM+1;
        %FM algorithm
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
    else
        countNW=countNW+1;
        if countNW ==1
            beginNW=countFM;
        end
        uIndex=0;
        [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
        xkArr=[xkArr;[xk',fk,1]];
    end
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
    if maxIter < countFM
        break;
    end
end
end
tf=etime(clock,t);
vk=sum(sign(rk));
disp(['%hybrid1 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);
%disp(['$',num2str(m),'\times ',num2str(n),'$&FM&(',num2str(countFM),',',num2str(countNW),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
%disp(['well1033&Daxs&',num2str(vk),'&',num2str(rn),'&',num2str(Ar),'&(',num2str(countFM),',',num2str(countNW),')&',num2str(beginNW)]);



function xk=MPRGP(A,b,x0,L,a,delta,maxIter)
k=0;r=A*x0-b;
tol=1e-15;
[m,n]=size(r);
[fx,bx]=fbxf(x0,m,tol,r);
p=fx;
bxn=bx'*bx;
while bxn+fx'*fx>delta
    f1x=(x0-l)./a;
    fxc=find(f1x>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx
       sqrt()
       disp(['conject gradient:',num2str(r'*r)]);
    else
       FM() 
       disp(['go to implement space:',num2str(r'*r)]);
    end
    k=k+1;
    bxn=bx'*bx;
    x0=xk;
end
end

function [amax,xc,res]=computeMaxStep(x0,d)
grt0=(d>0);
% if all indredient great 0 then natural sartisfied
if sum(grt0)<1
    amax=100;
    res=0;
else
    amax=min(x0(grt0)./d(grt0));
    res=1;
end
xc=x0+amax*d;
end

function [f,b]=fbxf(x0,m,tol,p)
AA=find(x0<tol);
FF=setdiff(1:m,AA);
b=zeros(m,1);
f=zeros(m,1);
b(AA)=p(AA);
f(FF)=p(FF);
b(b>0)=0;
end


function f=fxf(x0,m,tol,p)
FF=find(x0>tol);
f=zeros(m,1);
f(FF)=p(FF);
end

function b=bxf(x0,m,tol,p)
AA=find(x0<tol);
b=zeros(m,1);
b(AA)=p(AA);
b(b>0)=0;
end
