clc
clear
A=[1;-1];
b=[1;1];

x=5;
z=[1;1];
fm=1;

r=A*x-b-z;
f=0.5*r'*r;
g=[A'*r;-r];
p=[x;z]-g;
if p(2)<0
    p(2)=0;
end
if p(3)<0
    p(3)=0;
end

r1=A*p(1)-b-p([2;3]);
f1=0.5*r1'*r1;

AA=[A -eye(2)];
AAA=AA'*AA;
bbb=b'*AA;
0.5*p'*AAA*p-bbb*p+1
% A=[1,1;-1,-1;-1,0;-6,-3];
% b=[1;1;0.5;2];
% x0=[-1;-0.2];
% 
% 
% r0=A*x0-b;
% z0=r0;
% y0=-r0;
% AA=find(z0<1e-10);
% FF=find(z0>=1e-10);
% z0(AA)=0;
% y0(FF)=0;
% dF=A\y0;
% dN=A(AA,:)\(b(AA)-A(AA,:)*x0);
% AAAd=A(FF,:)*dN;
% zdN=-z0(FF)./AAAd;
% zdN(zdN<0)=inf;
% a1=min(zdN);
% zk=z0;
% xk=x0+a1*dN;
% zk(FF)=z0(FF)+a1*AAAd;
% 
% 
% AAAd=A*dN;
% zdN=-z0./AAAd;
% zdN(zdN<1e-10)=inf;
% a2=min(zdN);
% xk=x0+a2*dN;
% zk=z0+a2*AAAd;