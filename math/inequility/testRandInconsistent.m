m = 1000;
ratio=0.2;
n = floor( ratio * m);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
maxIter = 50;
nf = 3;
str = ['D','C','P','R'];
[xk,rk,countFM,countNW,beginNW,tf,vk,Arr,rkrr]=als(x0,A,b,maxIter);
fprintf('estiate value:%d\n',Arr(1,end));
iterA=size(4,maxIter+2);
maxIterA = 0;
for i=1:4
    type = str(i);
    [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,maxIter,nf,[type,'HA']);
    resvec = arvec;
    rkD=b-A*xkD;
    rkD(rkD<0)=0;
    dD=norm(rkD);
    gD=norm(A'*rkD);
    beginN=find(itersm>0);
    sumiter=sum(itersm>0);
    if isempty(beginN)
        beginN=0;
    end
    fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
    % fprintf('$ %d \\times %d $ & %g & %g & %4.4f & %g & %g & %g\n',m,n,dD,gD,tfD,iter*nf,sumiter,beginN);
    iterA(i,1:iter+1)=resvec;
   if maxIterA < iter+1
        maxIterA = iter+1;
    end
end
type=['r','g','k','c'];
beginp = 1;
figure
semilogy(beginp:maxIterA,iterA(1,beginp:maxIterA),'b.');
hold on
for i=2:4
    semilogy(beginp:maxIterA,iterA(i,beginp:maxIterA),[type(i) '.']);
end
legend('DHA','CHA','PHA','RHA');