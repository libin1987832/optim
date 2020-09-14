function [xs,rpk]=sm(A,b,n,rpk,x0)
AA=(rpk>0);
AI=A(AA,:);
bI=rpk(AA);
[u,flag,relres,iter,resvec,lsvec] = lsqrm(AI,bI,1e-12,n,[],[],zeros(n,1),A,b,x0,AA);
if flag == 5 || flag==1
    xs = x0 + u;
else
    aa = spiecewise(A,b,u,x0);
    xs = x0 + aa * u;
end