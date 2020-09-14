function [xs,rpk]=sm(A,b,tol,maxit,M1,M2,x0,ALL,bLL,xs,AA,varargin)
[u,flag,relres,iter,resvec,lsvec] = lsqrm(AI,bI,1e-12,n,[],[],zeros(n,1),A,b,x0,AA);
if flag == 5 || flag==1
    xs = x0 + u;
else
    aa = spiecewise(A,b,u,x0);
    xs = x0 + aa * u;
end