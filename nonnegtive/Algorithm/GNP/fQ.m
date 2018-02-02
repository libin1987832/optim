%罚函数方法的函数值
% fQ=1/2*||(Ax-b)+||^2+M||(-x)+|| i
function q=fQ(A,b,x0,M)
Y=A*x0-b;
Y(Y<0)=0;
XX=-1.*x0;
XX(XX<0)=0;
q=0.5*Y'*Y+M.*XX'*XX;
end
