function [x,err,index,indexN]=PA(x0,nmax,M,q)
tol=1e-10;
nf=5;
err=test_valid(M,q,x0);
index=1;
indexN=0;
while err>tol && index< nmax
    index=index+1;
    % Cauthy step
    [y,res]=fpi(x0,5,M,q);
    [a,xpf,result]=projectedsearch(y-x0,x0,M,q);
    I=(xpf>0);
    MII=M(I,I);
    qI=q(I);
    t=predict(MII,xpf(I),qI,4);
    t(t>0)=1;
    t(t<0)=0;
    cs=sum(I);
    if sum(t)>cs*0.985
        % subspace step
        xs=subspacesearch(xpf,M,q);
         [a,xk]=projectedsearch(xs-xpf,xpf,M,q);
        x0=xk;
        indexN=indexN+1;
    else
        x0=xpf;
    end
    % check for optimality
    err=test_valid(M,q,x0);
end
x=x0;