A=rand(6,5);
x0=rand(5,1);
XX=zeros(5,6);
YY=zeros(6,6);
XX(:,1)=x0;
for i=2:6
    YY(:,i-1)=A*XX(:,i-1);
    XX(:,i)=A'*YY(:,i-1);
end
C=XX(:,1:5)\XX(:,6);
C1=[[0;1;0;0;0] [0;0;1;0;0] [0;0;0;1;0] [0;0;0;0;1] C];
[qx rx]=qr(XX(:,1:5));
[qy ry]=qr(YY(:,1:5));
qy=qy(:,1:5);
ry=ry(1:5,:);
crx=C1'*rx'
iry=inv(ry')
% crx
% A=[1 2 3;2 4 7;3 7 5];
% x=[1;2;3];
% x1=A*x;
% x2=A*x1;
% x3=A*x2;
% XX=[x,x1,x2];
% C=XX\x3;
% C1=[[0;1;0] [0;0;1] C];
% [Q,R]=qr(XX);
% R*C1*inv(R)
