m = 10000;
ratio=0.1;
n = floor( ratio * m);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
maxIter = 10;
nf = 3;
str = ['D','C','P','R'];
[xk,rk,countFM,countNW,beginNW,tf,vk,Arr,rkrr]=als(x0,A,b,maxIter);
iterA=size(4,maxIter+2);
maxIterA = 0;
for i=1:4
    type = str(i);
    [xkD,flag,relres,iter,resvec,itersm,tfD]=hybridA(A,b,x0,maxIter,nf,[type,'HA']);
    rkD=b-A*xkD;
    rkD(rkD<0)=0;
    dD=norm(rkD);
    gD=norm(A'*rkD);
    beginN=find(itersm,true);
    sumiter=sum(itersm);
    fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN);
   % fprintf('$ %d \\times %d $ & %g & %g & %4.4f & %g & %g & %g\n',m,n,dD,gD,tfD,iter*nf,sumiter,beginN);
    iterA(i,1:iter+1)=resvec;
   if maxIterA < iter+1
        maxIterA = iter+1;
    end
end
type=['r','g','k','c'];
figure
semilogy(1:maxIterA,iterA(1,1:maxIterA),'b.');
hold on
for i=2:4
    semilogy(1:maxIterA,iterA(i,1:maxIterA),[type(i) '.']);
end
legend('DHA','CHA','PHA','RHA');