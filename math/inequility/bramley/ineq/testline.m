A = [-1;1;-1;-1];
b = [-4;7;-2;-7];
x = 0;
d = 1;
[stepsize,newfval,newderiv,err] = linesrch(A,x,d,r,fval,epsline,linemeth,neq,printlevel)