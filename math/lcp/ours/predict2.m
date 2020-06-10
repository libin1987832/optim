function [x]=predict2(x0,x1,x2)
d1=x1-x0;
d2=x2-x1;
lambda=max(abs(d2./d1));
lambda1=norm(d2)/norm(d1);
x=x0+(x1-x0)/(1-lambda);