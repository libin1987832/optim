m = 1000;
ratio=0.1;
n = floor( ratio * m);
A = 2 * rand(m , n)-1;
b = 2 * rand(m , 1)-1;
x0 = zeros(n , 1);
maxIter = 100;
nf = 3;
str = ['DHA''CHA','PHA','RHA'];
for i=1:4
    type = str(i);
    [xkD,flag,relres,iter,resvec,itersm,tf]=hybridA(A,b,x0,maxIter,nf,[type,'HA']);
    rkD=b-A*xkD;
    rkD(rkD<0)=0;
    dD=norm(rkD);
    gD=norm(A'*rkD);
    fprintf('%s $ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',type,m,n,dD,gD,tfD,countFD*nf,sum(countND,true),find(countND,true));
end