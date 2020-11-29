%计算非负函数的梯�?
% f=1/2*norm(Ax-b)+
% fdetq=A'(Ax-b)+
function g=fdetq(A,b,x0)
Y=b-A*x0;
Y(Y<0)=0;
g=-A'*Y;
end

