function [xk,rk]=Lei(x0,A,b,k,rk)
    rk(rk<0)=0;
    uk=krylovk(A,rk,k);
    xk=x0+uk;
    rk = b - A * xk;