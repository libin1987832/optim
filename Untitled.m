  n = 100; on = ones(n,1); A = spdiags([-2*on 4*on -on],-1:1,n,n);
       b = sum(A,2); tol = 1e-8; maxit = 15;
       M1 = spdiags([on/(-2) on],-1:0,n,n);
       M2 = spdiags([4*on -on],0:1,n,n);
       x = lsqr(A,b,tol,maxit,M1,M2);