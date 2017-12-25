function  a=fsearcha(A,b,x0,p)
t=0.5;
a0=1;
c1=0.01;
c2=0.9;
q1=fq(A,b,x0+a0*p);
q0=fq(A,b,x0);

g1=fdetq(A,b,x0+a0*p);
g0=fdetq(A,b,x0);

%while x0+a0*p<0||q1>=q0+c1*a0*g0'*p||g1'*p<=c2*g0'*p;
while q1>=q0+c1*a0*g0'*p||g1'*p<=c2*g0'*p;
   
    t=t^2;
    a0=a0*t;
    q1=fq(A,b,x0+a0*p);
 
    g1=fdetq(A,b,x0+a0*p);
end
a=a0;
end

