function  a=fwolfepowersearcha(A,b,x0,p)
%wolfepowerÏßĞÍËÑË÷
k0=0.5;
k1=0.5;
a0=k0;
c1=0.01;
c2=0.9;
q1=fq(A,b,x0+a0*p);
q0=fq(A,b,x0);

g1=fdetq(A,b,x0+a0*p);
g0=fdetq(A,b,x0);
while q1>=q0+c1*a0*g0'*p
    a0=k0*a0;
end
while g1'*p<=c2*g0'*p 
    b=(1/k0)*a0;
    a1=b;
    j=1;
    q3=fq(A,b,x0+a1*p);
    while q3>=q0+c1*a1*g0'*p
        a1=a0+k1^j*(b-a0);
        j=j+1;
    end
    a0=a1;
end
a=a0;
end

