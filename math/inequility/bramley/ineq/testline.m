A = [1;-1;1;1]*5;
b = [-4;3;-2;-7];
x0 = 0;
d = -1;
r = -b;
fval = 0.5*[0,3,0,0]*[0,3,0,0]';
linemeth = 1;
epsline = 1e-10;
neq = 0;
printlevel = 1;

x=[-5:0.1:5]';
n=size(x,1);
bM=repmat(b,1,n);
AM=repmat(A,1,n);
xM=diag(x);
rM=bM-AM*xM;
rMP = max(rM,zeros(size(rM)));
dd = 0.5*rMP'*rMP;
yy=diag(dd);
plot(x,yy);
 [stepsize,newfval,newderiv,err] = linesrch(-A,x0,d,r,fval,epsline,linemeth,neq,printlevel);
 addpath('../../.');
addpath('../../util');
[xk,flag,relres,iter,resvec,arvec,itersm,tf] = hybridA(A,b,x0,10,3,'DHA');
alph=spiecewise(A,b,-1,x0);
[stepsize alph,xk]
