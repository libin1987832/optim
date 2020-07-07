% syms x y
% %y=((1+0.1)^x)/(1+0.1*x);
% limit(((1+0.1)^x)/(1+0.1*x),x,inf)

clc
clear
addpath('../FM')
tfA1=[];
dfA1=[];
tfA6=[];
dfA6=[];
dim={};
maIter=500;
nmax=500;
etc=0.5;
ete=2;
rou=0.99;
trmax=1e12;
trr=1;

for m=1000:1000:2000
for ratio=0.1:0.01:0.3
        n=ceil(ratio*m);
        A=2*rand(m,n)-1;
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
        t=clock;
        [xk1,fk1,xkArr1,countF1,countN1]=hybrid1(x0,A,b,maIter);
