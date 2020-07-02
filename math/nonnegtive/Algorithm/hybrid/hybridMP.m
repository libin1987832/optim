function [xk,fk,xkArr,countFM,countNW,Q]=hybridMP(x0,A,b,maxIter)
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

%b1=zeros(n,1);
f1=zeros(n,1);
b2=zeros(m,1);
f2=zeros(m,1);
zk=zeros(m,1);
if Ar<delt*rn || rn<delt
    xk=x0;
    rk=r;
    fk=0.5*(r'*r);
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    bxn=bx'*bx;
    f1x=x0./a;
    fxc=find(f1x>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx
        countNW=countNW+1;
        if countNW ==1
            beginNW=countFM;
        end
        uIndex=0;
        %[xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
        z0=A*x0-b;
        AA=find(z0<ee);
        FF=setdiff(1:m,AA);
        d=-A(AA,:)\(z0(AA));
        %a1=-min(z./A*d);
        Ad=A(FF,:)*d;
        a2=-min(z(FF)./Ad);
        a2=min(a2,1);
        xk=x0+a2*d;
        zk=z0;
        zk(FF)=z0(FF)-a2*Ad;
        %zk=A*xk-b;
        %zk(zk<0)=0;
        Axkz=( A*xk-b-zk);
        f1=A'*Axkz;
        f2=0;
        f2(FF)=-Axkz(FF);
        b2=0;
        b2=-Axkz(AA);
        b2=min(b2,0);
        xkArr=[xkArr;[xk',fk,1]];
        disp(['conject gradient:',num2str(rk'*rk)]);
    else
        countFM=countFM+1;
        %FM algorithm
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        zk=-fm;
        zk(zk<0)=0;
        AA=find(zk<ee);
        FF=setdiff(1:m,AA);
        Axkz=(A*xk-b-zk);
        f1=A'*Axkz;
        f2=0;b2=0;
        f2(FF)=-Axkz(FF);
        b2=-Axkz(AA);
        b2=min(b2,0);
        xkArr=[xkArr;[xk',fk,0]];
        disp(['go to implement space:',num2str(rk'*rk)]);
    end
    bx=b2;
    fx=[f1;f2];
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



% function [f,b]=fbxf(x0,m,tol,p)
% AA=find(x0<tol);
% FF=setdiff(1:m,AA);
% b=zeros(m,1);
% f=zeros(m,1);
% b(AA)=p(AA);
% f(FF)=p(FF);
% b(b>0)=0;
% end
