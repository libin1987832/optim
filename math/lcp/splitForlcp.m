%function [x,err]=splitForlcp(x0,nmax,jc,je,delt0,deltmax,M,q)
function [x,err]=splitForlcp(x0,nmax,M,q)
tol=1e-5;
U=triu(M,1)
B=M-U;
nf=5;
err=test_valid(M,q,x0);
index=1;
while err>tol && index< nmax
 index=index+1;
 % Cauthy step
 [y,res]=fpi(x0,1,B,U,q);
 if res <0
     xk=y;
     break;
 end
 [ac,xc]=computeCauchyStep(q,M,x0,y-x0);
 % additional fixed point iterations
[y,res]=fpi(xc,nf,B,U,q);
 if res <0
     xk=y;
     break;
 end
 % projected search on fixed point iteration
 [a,xpf]=projectedsearch(y-xc,xc,M,q);
 % subspace step
 xs=subspacesearch(xpf,M,q);
 % projected search on subspace direction
 [a,xk]=projectedsearch(xs-xpf,xpf,M,q);
 x0=xk;
 % check for optimality
 err=test_valid(M,q,x0);
end
x=xk;