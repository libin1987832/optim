%function [x,err]=splitForlcp(x0,nmax,jc,je,delt0,deltmax,M,q)
function [x,err,index]=splitForlcp(x0,nmax,nf,M,q)
tol=1e-10;
err=test_valid(M,q,x0);
index=1;
while err>tol && index< nmax
 index=index+1;
 % Cauthy step
 [y,res]=fpi(x0,1,M,q);
 [ac,xc,result]=computeCauchyStep(q,M,x0,y-x0);
 % additional fixed point iterations
[y,res]=fpi(xc,nf,M,q);
 % projected search on fixed point iteration
 [a,xpf,result]=projectedsearch(y-xc,xc,M,q);
 % subspace step
 xs=subspacesearch(xpf,M,q);
 % projected search on subspace direction
 [a,xk]=projectedsearch(xs-xpf,xpf,M,q);
 x0=xk;
 % check for optimality
 err=test_valid(M,q,x0);
disp(['sym err:',num2str(err)]);
end
x=xk;