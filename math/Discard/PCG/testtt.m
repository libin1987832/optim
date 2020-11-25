%  A=[1 2 -1 0 0;2 1 0 -1  0;-1 -1 0 0 -1]
%   Ah=[1,2;2,1;-1,-1];
%   bk=[2;2;-1];
% [x1,f1,residual,exitflag]=lsqlin(A,bk,[],[],[],[],zeros(5,1),Inf*ones(5,1))
% y=A*x1-bk
% y1=A*(x1-[0.01,0.01,0,0,0]')-bk
% y2=A*[2/3;2/3;0;0;0]-bk
% norm(y)
% norm(y1)
% norm(y2)
% norm(bk)
% (Ah'*Ah)\(Ah'*bk)
% Ah'*(-bk)

N=2000;
p=rand([3,N])*100;
p=[p;ones(1,N)];
% p=[ 43.8744   79.5200   44.5586   75.4687   65.5098   49.8364   58.5268   25.5095   89.0903   13.8624
%    38.1558   18.6873   64.6313   27.6025   16.2612   95.9744   22.3812   50.5957   95.9291   14.9294
%    76.5517   48.9764   70.9365   67.9703   11.8998   34.0386   75.1267   69.9077   54.7216   25.7508
%    1          1          1         1         1       1            1        1        1         1    ];

yaw = 0.5854; 
pitch = 0.3; 
roll = 0.2;
dcm = angle2dcm( yaw, pitch, roll );
tran=[1,2,4];
T=[dcm,tran';0,0,0,1];
randy=0.01;
randp=0.01;
randr=0.01;
dcm0=angle2dcm( yaw+randy, pitch+randp, roll+randr );
T0=[dcm0,tran'+[0.01,0.02,0.01]'*0.01;0,0,0,1];
q=T*p;
q0=q(1:3,:)+rand(3,N)/100;
q0=[q0;ones(1,N)];
%J=[q(4,1)*eye(3,3),-1*[0 -q(3,1) q(2,1);q(3,1) 0 -q(1,1);-q(2,1) q(1,1) 0]];
J=[];
f=[];
tq=T0*p;
for i=1:N
    
    J=[J ; [eye(3,3),-1*[0 -tq(3,i) tq(2,i);tq(3,i) 0 -tq(1,i);-tq(2,i) tq(1,i) 0];]];
    pp=T0*p(:,i)-q0(:,i);
    f=[f;pp(1);pp(2);pp(3)];
end
detX=-(J'*J)\(J'*f);
detX=detX;
rot3f=detX(4:6);
tran3f=detX(1:3);
theta=norm(rot3f)
rot3fN=rot3f/theta;
itheta=1/theta;
A=sin(theta)*itheta;
B=(1-cos(theta))*(itheta*itheta);
C=(1-A)*(itheta*itheta);
cr=cross(rot3f,tran3f);
w_cr=cross(rot3f,cr);
trandetX=tran3f+B*cr+C*w_cr;
 rotdetX2=rotationVectorToMatrix(rot3f);
wx=rot3f(1);
wy=rot3f(2);
wz=rot3f(3);
wx2=wx^2;
wy2=wy^2;
wz2=wz^2;
rotdetX=[1-B*(wy2+wz2),B*(wx*wy)-A*wz,B*(wx*wz)+A*wy;
         B*(wx*wy)+A*wz,1-B*(wx2+wz2),B*(wy*wz)-A*wx;
         B*(wx*wz)-A*wy,B*(wy*wz)+A*wx,1-B*(wx2+wy2)];
T2=[rotdetX,trandetX;0 0 0 1]*T0;
[yaw2 pitch2 roll2]=dcm2angle( rotdetX*dcm0 );
% [yaw pitch roll ; yaw+randy, pitch+randp, roll+randr ;yaw2 pitch2 roll2]*180
norm(T*p-q0)
norm(T0*p-q0)
norm(T2*p-q0)

