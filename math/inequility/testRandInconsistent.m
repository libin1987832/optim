% for m = 1000:1000:2000
% for ratio = 0.1:0.1:0.3
    m = 1000;
    ratio=0.1;
    n = floor( ratio * m);
    A = 2 * rand(m , n)-1;
    b = 2 * rand(m , 1)-1;
    x0 = zeros(n , 1);
    maxIter = 50;
    nf = 3;
    str = ['D','C','R','P'];
    [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
     rkh=b-A*xkh;
     rkh(rkh<0)=0;
    dh=norm(rkh);
     gh=norm(A'*rkh);
    fprintf('han:%g %4.4f\n',gh,tfh);
%    [xk,rk,countFM,countNW,beginNW,tf,vk,Arr,rkrr]=als(x0,A,b,maxIter);
  %  fprintf('estiate value:%d\n',Arr(1,end));
    iterA=size(4,maxIter+2);
    maxIterA = 0;
    fprintA=zeros(1,8);
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
        fprintA(2*i-1:2*i)=[gD,tfD];
        if isempty(beginN)
            beginN=0;
        end
        fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
        iterA(i,1:iter+1)=resvec;
        if maxIterA < iter+1
            maxIterA = iter+1;
        end
    end
    fprintf('$ %d \\times %d $ & %g & %4.4f & %g & %4.4f & %g & %4.4f & %g & %4.4f\\\\\n'...
        ,m,n,fprintA(1),fprintA(2),fprintA(3),fprintA(4),fprintA(5),fprintA(6),fprintA(7),fprintA(8));
% end
% end

type=['r','g','k','c'];
typet=['+','o','v','s'];
beginp = 1;
figure
%maxIterA = 70;
h=semilogy(beginp:maxIterA,iterA(1,beginp:maxIterA),'bx');
h.LineStyle = '--';
hold on
for i=2:4
    h=semilogy(beginp:maxIterA,iterA(i,beginp:maxIterA),[type(i) typet(i)]);
    h.LineStyle = '--';
end
legend('DHA','CHA','RHA','PHA');
xlabel('Iteration Number');
ylabel('the norm of the gradient');
    