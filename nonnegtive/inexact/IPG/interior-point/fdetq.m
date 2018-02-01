%计算非负函数的梯度
% f=1/2*norm(Ax-b)+
% fdetq=A'(Ax-b)+
function g=fdetq(A,b,x0)
Y=A*x0-b;
Y(Y<0)=0;
g=A'*Y;
end

