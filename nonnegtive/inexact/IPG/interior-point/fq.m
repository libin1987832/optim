%非负最小二乘函数
% fq=||(Ax-b)+||
function q=fq(A,b,x0)
Y=A*x0-b;
Y(Y<0)=0;
q=0.5*(Y'*Y);
end

