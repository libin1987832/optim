function [y,res]=fpi(x0,nf,M,q)
[m,n]=size(M);
max_iter = 10;
tol_rel  = 0.0;
tol_abs  = 0.0;
for i=1:nf
%     [xk,s,iter,Aopt]=qp_bnd(B,d);
%    [w,z,retcode] = LCPSolve(B,d);
%     [z,res] = splitD(B,d,10);
%     [z,res] = splitS(M,q,1,x0,10);
    [z err iter flag convergence msg] =  pgs(M,q, x0, 1, tol_rel, tol_abs, false);
%     [z err iter flag convergence msg] = psor(M,q, x0, 1.4, 1, tol_rel, tol_abs, false)
%     x = LCP(B,d);
%    test_bnd(B,d,z)
    res=test_valid(M,q,x0);
    res1=test_valid(M,q,z);
    if res>res1
        x0=z;
    else
        res=-1;
        break;
    end
%     test_bnd(B,d,xk)
end
y=x0;


