function [xk,rk]=FMD(x0,A,b,rk,type)
rk(rk<0) = 0;
% compute min increase
if type ==1
    uk=A\rk;
else
    uk = lsqminnorm(A,rk);
end
xk=x0+uk;
rk=b-A*xk;