%function [x,err]=splitForlcp(x0,nmax,jc,je,delt0,deltmax,M,q)
function [x,err,index]=hybridSplit(x0,A,b,nmax,nf,M,q)
%function [x0,iter,nss]=asysub(x0,M,q,nmax,etc,ete,trr,trmax,rou)
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
nf=2;
ns=2;
tr0=(trr+trmax)/2;
iter=0;
nss=0;
while Ar>delt*rn && rn>delt
    cmax=0;
    iter=iter+nf+ns;
   [xfA,cmax] = splitS_asy_FM(A,b,1,x0,nf,cmax);
    
   [xs,zk]=kyrlov(AALL,b,x0);
   
    xsn=norm(xs-xfA(:,nf));
    mtr=min(1,tr0/xsn);
    xi=xfA(:,nf)+mtr*(xs-xfA(:,nf));
    xi=max(0,xi);
    nss=nss+ns;
    [xsA,cmax] = splitS(M,q,1,xi,ns,cmax);
    rok=max(rou,(1+cmax)/2);
    p1=norm(xsA(:,1)-xfA(:,nf));
    p2=norm(xsA(:,2)-xsA(:,1));
    p3=norm(xfA(:,nf)-xs(:,1));
    funfi3=funfi(xsA(:,ns),M,q);
    if p1<=rok*p3 && p2<=rok*p1
        x0=xsA(:,ns);
        tr0=median([trr,ete*tr0,trmax]);
    elseif funfi3<=0.5*funfim0
          x0=xsA(:,ns);
          funfim0=0.5*funfi3;
          tr0=median([trr,ete*tr0,trmax]);
    else
        tr0=etc*tr0;
    end
    if funfi3<1e-5
        break;
    end
end

tf=etime(clock,t);
vk=sum(sign(rk));
disp(['%hybridSplit m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);

end

function [xk,zk]=kyrlov(AALL,b,x0)
[m,n]=size(AALL);
fm0=AALL*x0-b;
z0=fm0;
z0(z0<0)=0;
ee=1e-15;% computer floating point arithmetic
AA=(z0<ee);
%FF=setdiff(1:m,AA)';

A=AALL(AA,:);

y=b(AA)-A*x0;

rmrk=AALL*x0-b-z0;
watch1=0.5*rmrk'*rmrk;


u1=0;  
beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;

ro_1=alph1;
thgma_1=beta1;
g1=v1;

for i=1:n
    q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
    
    ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
    
    v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
    
    theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
    
    d=thgma1*g1./ro1;
    u2=u1+d;  g2=v2-theta2*g1./ro1;
    
    xk=x0+u2;
    
    fmk=AALL*xk-b;
    zk=fmk;
    zk(zk<0)=0;
    rmrk=AALL*xk-b-zk;
    watch2=0.5*rmrk'*rmrk;

    AAk=(zk<ee);
    empty=isempty(setdiff(AA,AAk));
    if ~empty
        x1=x0+u1;
        fmk=AALL*x1-b;
        zk=fmk;
        zk(zk<0)=0;
        FF=setdiff(1:m,AA)';
        Ad=A(FF,:)*d;
        z0Ad=zk(FF)./Ad;
        z0Ad(z0Ad<=0)=inf;
        a2=min(z0Ad);
        a2=min(a2,1);
        xk=x1+a2*d;
        zk(FF)=zk(FF)-a2*Ad;
        
        rmrk=AALL*xk-b-zk;
        watch3=0.5*rmrk'*rmrk;

        break;
    end
    
    u1=u2;
    q1=q2;v1=v2;alph1=alph2;
    
    ro_1=ro_2;
    thgma_1=thgma_2;
    g1=g2;
end
end