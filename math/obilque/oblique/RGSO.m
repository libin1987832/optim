function [x,k]=RGSO(A,b,n,RRE,x,b_A)
%b_A=b-(A'\(A'*b));
RRE=RRE*norm(b,2)^(2);
N=zeros(n,1);%用于存储矩阵列范数
for i=1:1:n
    N(i)=norm(A(:,i),2)^(2);
end
%第一次迭代
i_1=ceil(rand(1,1)*n);
r=b-A*x;
alpha=A(:,i_1)'*r/N(i_1);
x(i_1)=x(i_1)+alpha;
r=r-alpha*A(:,i_1);
i_2=i_1;
for k=1:1:499999
    while i_1==i_2
        i_2=ceil(rand(1,1)*n);
    end
    G=A(:,i_1)'*A(:,i_2);
    g=N(i_2)-G^2/N(i_1);
    if g>0.0001
        alpha=A(:,i_2)'*r/g;
        beta=-G*alpha/N(i_1);
        x(i_2)=x(i_2)+alpha;
        x(i_1)=x(i_1)+beta;
        r=r-alpha*A(:,i_2)-beta*A(:,i_1);
    end
    if  norm(r-b_A)^2<RRE
        break;
    else
        i_1=i_2;
    end
end
k=k+1;
end
