m = 1000;
ratio=0.1;
n = floor( ratio * m);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
maxIter = 100;
nf = 3;
str = ['D','C','P','R'];
for i=1:4
    type = str(i);
    [xkD,flag,relres,iter,resvec,itersm,tfD]=hybridA(A,b,x0,maxIter,nf,[type,'HA']);
    rkD=b-A*xkD;
    rkD(rkD<0)=0;
    dD=norm(rkD);
    gD=norm(A'*rkD);
    beginN=find(itersm,true);
    sumiter=sum(itersm,true);
    fprintf('$ %d \\times %d $ & %g & %g & %4.4f & %g & %g & %g\n',m,n,dD,gD,tfD,iter*nf,sumiter,beginN);
end