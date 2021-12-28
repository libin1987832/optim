% A=qr(rand(3,3));
% [q,r]=qr(A)
% T=q'*diag([1/2,1/3,1])*q;
% T'*T
% A=sym([0.5,0.3;1,0.7;1,3;0.3,0.6])
% [q,r]=qr(A)
% Aq=q(:,[1,2]);
% B=diag([1,1,1,1])-Aq*Aq'*diag([1,1,0,0])
% [v,j]=jordan(B);
% vpa(j,4)
% tt=v*diag([0,0,1,1])*inv(v)
% vpa(tt,4)
% tt=v*inv(v)
% vpa(tt,4)
N=1000;
for i =1:1000
A=rand(N,N);
x=rand(N,1);
y=rand(N,1);
z=x+y;
z(z<0)=0;
f1=norm(A*z)^2;
x(x<0)=0;
y(y<0)=0;
f2=norm(A*x)^2+norm(A*y)^2;
if 2*f2-f1<0
   2*f2-f1
    break;
end
end
i