%罚函数方法的梯度
% f=1/2*||(Ax-b)+||^2+M||(-x)+|| 
% det1F=A'(Ax-b)+-2M(-x)+
function v=det1F(x,A,b,M)
    B=A*x-b;
    B(B<0)=0;%(Ax-b)+
    x0=-x;
    x0(x0<0)=0;%(-x)+
    v=A'*B-2*M*x0;
