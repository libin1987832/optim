function [xk,r0,rk,fk,fm,fr]=FMn(x0,Q,R,A,b,n)
for i=1:n
    %FM algorithm
    [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
    x0=xk;
end
