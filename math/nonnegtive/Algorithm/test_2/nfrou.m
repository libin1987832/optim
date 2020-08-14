clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')


for batch=1:10
Arecord=[];
%     bnf=50;
% enf=55;
% nnf=enf-bnf+1;
    for m=200:400:1000
    for ratio=0.1:0.1:1
%         for nf=bnf:enf
        n=floor(ratio*m);
         %A=2*rand(m,n)-1;

%          A = 2*randn(m,m) * ...
%          [diag([ones(floor(n/2),1)*100;ones(ceil(n/2),1)])*10;ones(m-n,n)*10] * ...
%          randn(n,n)-1;
     
         b=2*rand(m,1)-1;
%          load('test')
        x0=zeros(n,1);
        xr=rand(n,1);
        
        t1=etime(clock,t);
        t=clock;
        pinv(A);
        t2=etime(clock,t);
        fprintf('time A*x:%g,pinv:%g,qr:%g,r*q*x:%g\n',t1,t2,tq,tr);
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
         fprintf('ALS$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d,%g\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end),t1*countFA/10);
         Arecord=[Arecord;m n Arr(1,end) gA];
        
end
end
save(['ff' num2str(batch) '.mat'],'Arecord')
end
%%
clear
f1=load('ff1.mat');
f2=load('ff2.mat');
f3=load('ff3.mat');
f4=load('ff4.mat');
f5=load('ff5.mat');
f6=load('ff6.mat');
f7=load('ff7.mat');
f8=load('ff8.mat');
f9=load('ff9.mat');
f10=load('ff10.mat');
summ=[f1.Arecord(:,1:3),f2.Arecord(:,3),f3.Arecord(:,3),f4.Arecord(:,3),...
    f5.Arecord(:,3),f6.Arecord(:,3),f7.Arecord(:,3),f8.Arecord(:,3)...
    ,f9.Arecord(:,3),f10.Arecord(:,3)];
summ2=[summ mean(summ(:,3:12),2) std(summ(:,3:12),0,2)];
subplot(1,3,1)
plot(0.1:0.1:1,summ2(1:10,13))
subplot(1,3,2)
plot(0.1:0.1:1,summ2(11:20,13))
subplot(1,3,3)
plot(0.1:0.1:1,summ2(21:30,13))

