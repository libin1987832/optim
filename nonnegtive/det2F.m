%罚函数的二阶梯度
% f=1/2*||(Ax-b)+||^2+M||(-x)+|| i
% det2F=A'DA+2MD0
function v=det2F(x,A,b,M)
    B=A*x-b;
    D=diag(B);
    D(B>0)=1;
    D(B<0)=0;
    D0=diag(x);
    D0(D0>0)=0;
    D0(D0<0)=1;
    v=A'*D*A+2*M*D0;
