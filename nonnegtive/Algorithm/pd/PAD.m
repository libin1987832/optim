function [xk,fk]=PAD(x0,A,b)
y=A*x0-b;
All=1:n;
NE=All(y<=0);
Aa=A(NE,:);
Ba=b(NE,:);
[xk,fk]=kyrlov(x0,Aa,Ba);