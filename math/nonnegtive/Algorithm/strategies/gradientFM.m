function [xk,rk,countFM,countNW,beginNW,tf,vk]=gradientFM(x0,A,b,a,L,maxIter)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
%FM need a qr decompose
[Q,R]=qr(A);
rkp=b-A*x0;
r=rkp;
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
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
z0=-rkp;
z0(z0<0)=0;
AA=find(z0<ee);
%FF=setdiff(1:m,AA)';
        
Axkz=(-rkp-z0);
f1=A'*Axkz;
%f2(FF)=-Axkz(FF); will be zeros
% b2=-Axkz(AA);will be zeros
% b2=min(b2,0);will be zeros
bx=b2;
fx=[f1;f2];
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    bxn=bx'*bx;
    f1x=[x0;z0]./a;
    fxc=find(f1x>fx);
    f1x(fxc)=fx(fxc);
    if bxn<=L*f1x'*fx
        countNW=countNW+1;
        if countNW ==1
            beginNW=countFM;
        end
        uIndex=0;
        %[xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
        %d=-A(AA,:)\(fm0(AA));
        [xk,zk]=krylov(A,b,x0,rkp);
        Axkz=A*xk-b-zk;
        f1=A'*Axkz;
        f2(:)=0;
        AA=zk<ee;
        FF=setdiff(1:m,AA)';
        f2(FF)=-Axkz(FF);
        b2(:)=0;
        b2=-Axkz(AA);
        b2=min(b2,0);
        rk=b-A*xk;
        rk(rk<0)=0;
%         fk=0.5*rk'*rk;
%         xkArr=[xkArr;[xk',fk,1]];
 %       disp(['conject gradient:',num2str(fk)]);
    else
        countFM=countFM+1;
        %FM algorithm
        %[xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        [xk,r0,rk]=FixedM(x0,Q,R,A,b);
        zk=-rk;
        zk(zk<0)=0;
        AA=find(zk<ee);
        FF=setdiff(1:m,AA)';
        Axkz=(A*xk-b-zk);
        f1=A'*Axkz;
        f2(:)=0;b2(:)=0;
        b2=-Axkz(AA);
        b2=min(b2,0);
        rk(rk<0)=0;
        %xkArr=[xkArr;[xk',fk,0]];
%        disp(['go to implement space:',num2str(rk'*rk)]);
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

tf=etime(clock,t);
vk=sum(sign(rk));
