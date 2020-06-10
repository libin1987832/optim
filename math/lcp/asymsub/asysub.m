function [x0,iter,nss]=asysub(x0,M,q,nmax,etc,ete,trr,trmax,rou)
funfi0=funfi(x0,M,q);
funfim0=max(funfi0,1e5);
nf=2;
ns=2;
tr0=(trr+trmax)/2;
iter=0;
nss=0;
for i=1:nmax
    cmax=0;
    iter=iter+nf+ns;
   [xfA,cmax] = splitS_asy(M,q,1,x0,nf,cmax);
    xs=subspacesearch(xfA(:,nf),M,q);
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
x0;