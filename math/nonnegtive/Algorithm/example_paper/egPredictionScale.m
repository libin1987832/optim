addpath('../FM')
addpath('../hybrid')
statistic=0;
for m=100:50:1000
    ratio=0.1;
    n=ceil(ratio*m);
    A=2*rand(m,n)-1;
    b=2*rand(m,1)-1;
    x0=zeros(n,1);
    [xk1,fk1,xkArr1,countF1,countN1]=hybrid6(x0,A,b,20);
    if countN1==1
        statistic=statistic+1;
    end
end

