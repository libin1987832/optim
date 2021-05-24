function [nf,iffind]=generateNf(A,b,x0,iter)
rk0 = b-A*x0;
iter = 50;
rkA=zeros(m,iter);
sumrkA=zeros(1,iter);
rkA(:,1)=rk0;
sumrkA(1)=sum(rk0>0);
rk=rk0;
for i = 2:iter
[xk,rk]=Lei(x0,A,b,3,rk);
x0=xk;
rkA(:,i)=rk;
sumrkA(i)=sum(rk>0);
end
for i=1:iter
    signrks = rk>0;
    signi = rkA(:,i)>0;
    alls = all(signrk-signi);
end

