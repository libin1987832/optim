function test_fpi(x0,xk,B,C,q)
 d=C*x0+q;
 res=min(xk,B*xk+d);
 err=max(res);
 
 A=B+C;
f0=0.5*x0'*A*x0+q'*x0;
fk=0.5*xk'*A*xk+q'*xk;



% test example
