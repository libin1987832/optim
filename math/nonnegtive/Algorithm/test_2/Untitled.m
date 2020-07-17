% A=[1 0 0;2 1 0;3 0 1];
% B=[1 0 0;0 1 0;0 4 1];
% inv(A)*inv(B)
% A*B
% B*A 
A=[1;-1];
b=[-1;2];
[Q,R]=qr(A);
x0=1;
r=b-A*x0;
[xk,rk]=FixedM(x0,Q,R,A,b,r)
[xk,rk]=FixedM(xk,Q,R,A,b,rk)
[xk,rk]=FixedM(xk,Q,R,A,b,rk)
[xk,rk]=FixedM(xk,Q,R,A,b,rk)