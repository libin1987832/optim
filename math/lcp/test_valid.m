function resn=test_valid(Q,q,x0)
Qd=Q*x0+q;
res=min(x0,Qd);
resn=norm(res,'inf');
