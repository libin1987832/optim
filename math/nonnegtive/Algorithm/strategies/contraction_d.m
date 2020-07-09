%function [x,err]=splitForlcp(x0,nmax,jc,je,delt0,deltmax,M,q)
function [xk,rk,countFM,countNW,beginNW,tf,vk]=contraction_d(x0,A,b,maxIter,nf,ns,etc,ete,trr,trmax,rou)
%function [x0,iter,nss]=asysub(x0,M,q,nmax,etc,ete,trr,trmax,rou)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
%FM need a qr decompose
[Q,R]=qr(A);
rkp=b-A*x0;
r=rkp;
r(r<0)=0;
funfim0=0.5*(r'*r);
%condition for terminate
Ar=norm(A'*r);
rn=norm(r);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;

countFM=0;
countNW=0;
beginNW=nf;

if Ar<delt*rn || rn<delt
    xk=x0;
    rk=r;
    fk=0.5*(r'*r);
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
tr0=(trr+trmax)/2;
iter=0;
nss=0;
while Ar>delt*rn && rn>delt
    cmax=0;
    iter=iter+nf+ns+1;
    countFM=iter;
    [xfA,rkp,cmax] = splitS_asy_FM(A,b,Q,R,x0,nf+1,cmax,rkp);
    countNW=countNW+1;
    
    [xs,zk]=krylov(A,b,x0,rkp);
    rkp=b-A*xs;
    xsn=norm(xs-xfA(:,nf));
    mtr=min(1,tr0/xsn);
    xi=xfA(:,nf)+mtr*(xs-xfA(:,nf));
    nss=nss+ns;
    %[xsA,cmax] = splitS(A,q,xi,ns,cmax);
    [xsA,rkp,cmax] = splitS_asy_FM(A,b,Q,R,xi,ns,cmax,rkp);
    
    rok=max(rou,(1+cmax)/2);
    p1n=xsA(:,1)-xfA(:,nf);
    p1=p1n'*p1n;
    p2n=xsA(:,2)-xsA(:,1);
    p2=p2n'*p2n;
    p3n=xfA(:,nf-1)-xfA(:,nf);
    p3=p3n'*p3n;
    rsk=rkp;
    rsk(rsk<0)=0;
  %  funfi3=funfi(xsA(:,ns),M,q);
    funfi3=0.5*(rsk'*rsk);
    %%
    fprintf('p1/p3:%g,p2/p1:%g,ro:%g,f:%g,f0:%g,tr:%g\n',p1/p3,p2/p1,rok,funfi3,funfim0,tr0);
    %%
    if p1<=rok*p3 && p2<=rok*p1
        xk=xsA(:,ns);
        tr0=median([trr,ete*tr0,trmax]);
    elseif funfi3<=0.5*funfim0
        xk=xsA(:,ns);
        funfim0=0.5*funfi3;
        tr0=median([trr,ete*tr0,trmax]);
    else
        tr0=etc*tr0;
        xk=x0;
    end
    rk=rsk;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;  
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));
%disp(['%hybridSplit m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);
