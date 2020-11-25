clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')
tfA1=[];
dfA1=[];
tfA6=[];
dfA6=[];
dim={};
maxIter=500;
nmax=500;
etc=0.5;
ete=2;
rou=0.99;
trmax=1e2;
trr=1;


for batch=1:10
Arecord=[];
    for m=100:50:1000
        n=100;

         A = 2*randn(m,m) * ...
         [diag([ones(floor(n/2),1)*100;ones(ceil(n/2),1)])*10;ones(m-n,n)*10] * ...
         randn(n,n)-1;
     
         b=2*rand(m,1)-1;
%          load('test')
        x0=zeros(n,1);
        xr=rand(n,1);
        xs=-1;
        
        [xkR,xkR2,countFR,countNWR,bNWR,tfR,vkR]=residualR(x0,A,b,maxIter);
         dR=norm(xkR-xs);
         rkR=b-A*xkR;
         rkR(rkR<0)=0;
         gR=norm(A'*rkR);
         dR2=norm(xkR2-xs);
         rkR2=b-A*xkR2;
         rkR2(rkR2<0)=0;
         gR2=norm(A'*rkR2);
         fprintf('resdual$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,gR2,gR,tfR,countFR,countNWR);
         
         [xkA,rkA,countFA,countNA,bNWA,tfA,vkA,Arr]=als(x0,A,b,maxIter);
         dA=norm(xkA-xs);
         rkA=b-A*xkA;
         rkA(rkA<0)=0;
         gA=norm(A'*rkA);
         fprintf('ALS$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d,%g\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end),countFA/10);
         Arecord=[Arecord;m n Arr(1,end) gA];
end
save(['ffdimesion' num2str(batch) '.mat'],'Arecord')
end

%%
clear
f1=load('ffdimesion1.mat');
f2=load('ffdimesion2.mat');
f3=load('ffdimesion3.mat');
f4=load('ffdimesion4.mat');
f5=load('ffdimesion5.mat');
f6=load('ffdimesion6.mat');
f7=load('ffdimesion7.mat');
f8=load('ffdimesion8.mat');
f9=load('ffdimesion9.mat');
f10=load('ffdimesion10.mat');
summ=[f1.Arecord(:,1:3),f2.Arecord(:,3),f3.Arecord(:,3),f4.Arecord(:,3),...
    f5.Arecord(:,3),f6.Arecord(:,3),f7.Arecord(:,3),f8.Arecord(:,3)...
    ,f9.Arecord(:,3),f10.Arecord(:,3)];
summ2=[summ mean(summ(:,3:12),2) std(summ(:,3:12),0,2)];
plot(100:50:1000,summ2(:,13))