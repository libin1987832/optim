function [y,res]=fpi(x0,nf,B,C,q)
[m,n]=size(B);
for i=1:nf
    d=C*x0+q;
%     [xk,s,iter,Aopt]=qp_bnd(B,d);
%    [w,z,retcode] = LCPSolve(B,d);
%     [z,res] = splitD(B,d,10);
    [z,res] = splitS(B,d,1,10);
%     x = LCP(B,d);
%    test_bnd(B,d,z)
    res=test_valid(B,d,x0);
    res1=test_valid(B,d,z);
    if res>res1
        x0=z;
    else
        res=-1;
        break;
    end
%     test_bnd(B,d,xk)
end
y=x0;


