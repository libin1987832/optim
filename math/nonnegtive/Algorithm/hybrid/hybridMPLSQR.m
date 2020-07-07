function [xk,fk,xkArr,countFM,countNW,Q]=hybridMPLSQR(x0,A,b,a,L,maxIter)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
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
z0=A*x0-b;
z0(z0<0)=0;
AA=find(z0<ee);
%FF=setdiff(1:m,AA)';
        
Axkz=(A*x0-b-z0);
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
        [xk,zk]=kyrlov(A,b,x0);
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
        fk=0.5*rk'*rk;
        xkArr=[xkArr;[xk',fk,1]];
 %       disp(['conject gradient:',num2str(fk)]);
    else
        countFM=countFM+1;
        %FM algorithm
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        zk=-fm;
        zk(zk<0)=0;
        AA=find(zk<ee);
        FF=setdiff(1:m,AA)';
        Axkz=(A*xk-b-zk);
        f1=A'*Axkz;
        f2(:)=0;b2(:)=0;
        b2=-Axkz(AA);
        b2=min(b2,0);
        xkArr=[xkArr;[xk',fk,0]];
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
disp(['%hybridMPLSQR m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);
%disp(['$',num2str(m),'\times ',num2str(n),'$&FM&(',num2str(countFM),',',num2str(countNW),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
%disp(['well1033&Daxs&',num2str(vk),'&',num2str(rn),'&',num2str(Ar),'&(',num2str(countFM),',',num2str(countNW),')&',num2str(beginNW)]);


% function [xk,zk]=kyrlov(AALL,b,x0)
% [m,n]=size(AALL);
% fm0=AALL*x0-b;
% z0=fm0;
% z0(z0<0)=0;
% ee=1e-15;% computer floating point arithmetic
% AA=(z0<ee);
% %FF=setdiff(1:m,AA)';
% 
% A=AALL(AA,:);
% 
% y=b(AA)-A*x0;
% 
% rmrk=AALL*x0-b-z0;
% watch1=0.5*rmrk'*rmrk;
% 
% 
% u1=0;  
% beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;
% 
% ro_1=alph1;
% thgma_1=beta1;
% g1=v1;
% 
% for i=1:n
%     q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
%     
%     ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
%     
%     v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
%     
%     theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
%     
%     d=thgma1*g1./ro1;
%     u2=u1+d;  g2=v2-theta2*g1./ro1;
%     
%     xk=x0+u2;
%     
%     fmk=AALL*xk-b;
%     zk=fmk;
%     zk(zk<0)=0;
%     rmrk=AALL*xk-b-zk;
%     watch2=0.5*rmrk'*rmrk;
% 
%     AAk=(zk<ee);
%     empty=isempty(setdiff(AA,AAk));
%     if ~empty
%         x1=x0+u1;
%         fmk=AALL*x1-b;
%         zk=fmk;
%         zk(zk<0)=0;
%         FF=setdiff(1:m,AA)';
%         Ad=A(FF,:)*d;
%         z0Ad=zk(FF)./Ad;
%         z0Ad(z0Ad<=0)=inf;
%         a2=min(z0Ad);
%         a2=min(a2,1);
%         xk=x1+a2*d;
%         zk(FF)=zk(FF)-a2*Ad;
%         
%         rmrk=AALL*xk-b-zk;
%         watch3=0.5*rmrk'*rmrk;
% 
%         break;
%     end
%     
%     u1=u2;
%     q1=q2;v1=v2;alph1=alph2;
%     
%     ro_1=ro_2;
%     thgma_1=thgma_2;
%     g1=g2;
% end
% end