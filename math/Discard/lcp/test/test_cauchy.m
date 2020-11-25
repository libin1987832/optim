function dd=test_cauthy(a,d,A,c,x0)
dd=x0+a*d;
inv=0:0.01*(a+0.3):a+0.3;
fA=[];
for i=1:length(inv)
    xk=x0+inv(i)*d;
    f=0.5*xk'*A*xk+c'*xk;
    fA=[fA;f];
end
plot(inv,fA);
