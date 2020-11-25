function test_bnd(M,q,x)
Mq=M*x+q; 
res=min(x,Mq);
err=max(res);
