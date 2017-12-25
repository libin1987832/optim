function v=det1F(x,A,b,M)
    B=A*x-b;
    B(B<0)=0;
    x0=-x;
    x0(X0<0)=0;
    v=A'*B-2*M*x0;