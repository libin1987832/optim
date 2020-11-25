%function [x,err]=splitForlcp(x0,nmax,jc,je,delt0,deltmax,M,q)
function [x,err,index]=splitForlcp1(x0,nmax,nf,M,q,iDLU,iDL)
tol=1e-10;
err=test_valid(M,q,x0);
index=1;
while err>tol && index< nmax
    index=index+1;
    % Cauthy step
    if index <3
        [y,res]=fpi(x0,1,M,q);
    else
        y=Gauss_Siedel(iDLU,iDL,q,x0,1);
    end
    [ac,xc,result]=computeCauchyStep(q,M,x0,y-x0);
    % additional fixed point iterations
    if index <3
        [y,res]=fpi(xc,nf,M,q);
    else
        y=Gauss_Siedel(iDLU,iDL,q,x0,nf);
    end
    % projected search on fixed point iteration
    [a,xpf,result]=projectedsearch(y-xc,xc,M,q);
    % subspace step
    xs=subspacesearch(xpf,M,q);
    % projected search on subspace direction
    [a,xk]=projectedsearch(xs-xpf,xpf,M,q);
    x0=xk;
    % check for optimality
    err=test_valid(M,q,x0);
    % disp(['sym err:',num2str(err)]);
end
x=xk;