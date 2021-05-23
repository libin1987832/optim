function [xk,rk]=FMLSQR(x0,A,b,rk)
rk(rk<0)=0;
uk=krylovk(A,rk,3);
xk=x0+uk;
rk = b - A * xk;