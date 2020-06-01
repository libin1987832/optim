function x0=asysub(x0,M,q,nmax,etc,ete,trr,trmax,rou)
funfi0=funfi(x0,M,q);
funfim0=max(funfi0,1e5);
nf=5;
ns=3;
tr0=(trr+trmax)/2;
for i=1:nmax
    [xfn,cmax] = splitS(M,q,1,x0,nf-1,cmax);
    [xfn1,cmax] = splitS(M,q,1,xfn,1,cmax);
    xs=subspacesearch(xfk,M,q);
    xsn=norm(xs-xk);
    mtr=min(1,tr0/xsn);
    xi=xi+mtr*(xs-xk);
    [xsk1,cmax] = splitS(M,q,1,xi,1,cmax);
    [xsk2,cmax] = splitS(M,q,1,xsk1,1,cmax);
    [xsk3,cmax] = splitS(M,q,1,xsk2,ns-1,cmax);
    rok=max(rou,(1+cmax)/2);
    p1=norm(xsk1-xfn1);
    p2=norm(xsk2-xsk1);
    p3=norm(xfn-xfn1);
    funfi3=funfi(xk3,M,q);
    if p1<=rok*p3 && p2<=rok*p1
        x0=xsk3;
        tr0=median([trr,ete*tr0,trmax]);
    elseif funfi3<=0.5*funfim0
          x0=xsk3;
          funfim0=0.5*funfi3;
          tr0=median([trr,ete*tr0,trmax]);
    else
        tr0=etc*tr0;
    end
    if funfi3<1e-6
        break;
    end
end
x0;