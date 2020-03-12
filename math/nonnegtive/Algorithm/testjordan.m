% test jordan
qI1=[1;-1;1;6];
q1=qI1/sqrt(39);
qI2=[19;-19;-20;-3];
q2=qI2/sqrt(1131);
Q=[q1,q2];
Q1=qI1*qI1'*1131;
Q2=qI2*qI2'*39;
rT=diag([39*1131,39*1131,39*1131,39*1131])-(Q1+Q2)*diag([1,1,0,0]);
rts=sym(rT/(39*1131))
% rtsn=diag([1,1,0,0])
[v,j]=jordan(rts);
rtsj=v*j*inv(v)
r=[3/2;1/2;-1/4;-7/4];
rs=sym(r);
rst=v*j*inv(v)*rs