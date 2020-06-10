function [resn,fx]=test_valid(Q,q,x0)
Qd=Q*x0+q;
fx=[x0,Qd];
res=min(x0,Qd);
resn=norm(res,'inf');
