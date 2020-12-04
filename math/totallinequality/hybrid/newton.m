function [xs,rpk,len,flag] = newton(A,b,n,rpk,x0)
    AA = (rpk>tol);
        RR = (x0>0);
        % subspace
        AI = A(AA,RR);
        bI = rpk(AA);
        u = lsqminnorm(AI,bI);
        p = zeros(n,1);
        p(RR) = u;
 %       debug for test.m nf = 2
        [alpha, aranges, retcode] = arraySpiece(A,b,x0,p);
        x0 = x0 + alpha*p;
        indexsm = indexsm + 1;
    