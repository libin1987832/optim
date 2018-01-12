function q=fQ(A,b,x0,M)
Y=A*x0-b;
Y(Y<0)=0;
XX=-1.*x0;
XX(XX<0)=0;
q=0.5*norm(Y,2)^2+M.*norm(XX,2)^2;
end
