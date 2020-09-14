clear;
clc;
m = 100;
ratio=0.5;
n = floor( ratio * m);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
maxIter = 100;
nf = 3;
str = ['D','C','P','R'];
 [xk,rk,countFM,countNW,beginNW,tf,vk,Arr,rkrr]=als(x0,A,b,maxIter);
 Arr(:,end)
iterA=size(4,maxIter);
 for i=1:4
    type = str(i);
    [xkD,flag,relres,iter,resvec,itersm,tfD]=hybridA(A,b,x0,maxIter,nf,[type,'HA']);
    iterA(i,1:iter+1)=resvec;
    rkD=b-A*xkD;
    rkD(rkD<0)=0;
    dD=norm(rkD);
    gD=norm(A'*rkD);
    beginN=find(itersm,true);
    sumiter=sum(itersm);
    fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN);
   % fprintf('$ %d \\times %d $ & %g & %g & %4.4f & %g & %g & %g\n',m,n,dD,gD,tfD,iter*nf,sumiter,beginN);
 end
type=['r','g','d','k'];
figure
semilogy(1:maxIter,iterA(1,:),'b.');
hold on
for i=2:4
    semilogy(1:maxIter,iterA(i,:),[type(i) 'o']);
end
    