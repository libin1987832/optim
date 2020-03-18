% test jordan A=[1,1;-1,-1;-1,0;-6,-3] b=[1;1;1/2;2];
qI1=[1;-1;-1;-6];
q1=qI1/sqrt(39);
qI2=[19;-19;20;3];
q2=qI2/sqrt(1131);
Q=[q1,q2];
Q1=qI1*qI1'*1131;
Q2=qI2*qI2'*39;
rT=diag([39*1131,39*1131,39*1131,39*1131])-(Q1+Q2)*diag([1,1,0,0]);
rts=sym(rT/(39*1131))
% rtsn=diag([1,1,0,0])
[v,j]=jordan(rts);
rtsj=v*j*inv(v)
v*diag([0,1,1,1])*inv(v)
x0=[-3/4;1/4+1/4];
% x0=[-3/4;1/4];
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;1/2;2];
r=b-A*x0;
% r=[3/2;1/2;-1/4;-7/4];
rs=sym(r);
rst=v*j*inv(v)*rs
v*diag([0,1,1,1])*inv(v)
v*diag([0,1,1,1])*inv(v)*rs
v*diag([0,1,1,1])*inv(v)*rst
nn1=null(v*diag([0,1,1,1])*inv(v));
nn1=nn1/norm(nn1);
nn2=(rs-rst);
nn2=nn2/norm(nn2);
% vpa(nn1)
% vpa(nn2)
vpa(rst);