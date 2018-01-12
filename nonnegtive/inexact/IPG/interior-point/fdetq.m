function g=fdetq(A,b,x0)
Y=A*x0-b;
Y(Y<0)=0;
g=A'*Y;
end

