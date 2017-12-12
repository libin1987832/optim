function  a=fsearchaM(A,b,x0,p,M)
t=0.5;
a0=1;
c1=0.01;
c2=0.9;
q1=fQ(A,b,x0+a0*p,M);
q0=fQ(A,b,x0,M);

g1=det1F(x0+a0*p,A,b,M);
g0=det1F(x0,A,b,M);

%while x0+a0*p<0||q1>=q0+c1*a0*g0'*p||g1'*p<=c2*g0'*p;
while q1>=q0+c1*a0*g0'*p||g1'*p<=c2*g0'*p;
   
    t=t^2;
    a0=a0*t;
    q1=fq(A,b,x0+a0*p);
 
    g1=fdetq(A,b,x0+a0*p);
end
a=a0;
end

