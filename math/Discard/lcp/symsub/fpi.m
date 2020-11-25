function [y,err]=fpi(x0,nf,M,q)
[m,n]=size(M);

tol_rel  = 1e-5;
tol_abs  = 1e-10;

[y err iter flag convergence msg] =  pgs(M,q, x0, nf, tol_rel, tol_abs, true);


