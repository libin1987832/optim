function test_projected(a,A,c,x0,d)
inv=0:0.01*a:2*a;
fA=[];
fA1=[];
for i=1:length(inv)
    xk=x0+inv(i)*d;
    f=0.5*xk'*A*xk+c'*xk;
    fA1=[fA1;f];
    xk(xk<0)=0;
    f=0.5*xk'*A*xk+c'*xk;
    fA=[fA;f];
end

plot(inv,fA);
hold on 
plot(inv,fA1);