function [nf,iffind,xk,sumrkA]=generateNf(A,b,x0,iter,k)
[m,n]=size(A);
rk0 = b-A*x0;
%iter = 50;
rkA=zeros(m,iter);
sumrkA=zeros(1,iter);
rkA(:,1)=rk0;
sumrkA(1)=sum(rk0>0);
rk=rk0;
for i = 2:iter
[xk,rk]=Lei2(x0,A,b,k,rk);
x0=xk;
rkA(:,i)=rk;
sumrkA(i)=sum(rk>0);
end
for i=1:iter
    signrks = rk>0;
    signi = rkA(:,i)>0;
    alls = all(~abs(signrks-signi));
    if alls == 1
        nf = i;
        iffind = 1;
        break;
    end 
end

