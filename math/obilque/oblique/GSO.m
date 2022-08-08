function [x,error]=GSO(A,b,RRE,x0,tol,iter,b_A,debug)
[~,n]=size(A);
%b_A=b-(A'\(A'*b));
RRE=RRE*norm(b,2)^(2);
error=zeros(1,iter+1);
%N=zeros(n,1);%用于存储矩阵列范数
% for i=1:1:n
%     N(i)=norm(A(:,i),2)^(2);
% end
N=sum(A.*A,1);
%第一次迭代
r=b-A*x0;
if debug
    error(1)=r'*r;
end
x=x0;
alpha=A(:,1)'*r/N(1);
x(1)=x(1)+alpha;
r=r-alpha*A(:,1);
i_1=1;
for k=1:1:iter
    i_2=mod(k,n)+1;
    G=A(:,i_1)'*A(:,i_2);
    g=N(i_2)-G^2/N(i_1);
    if g>tol
        alpha=A(:,i_2)'*r/g;
        beta=-G*alpha/N(i_1);
        ab=[alpha;beta];
        x([i_2,i_1],1)=x([i_2,i_1],1)+ab;
        r=r-[A(:,i_2),A(:,i_1)]*ab;
        if debug
            error(k+1)=r'*r;
        end
    end
    if  norm(r-b_A)^2<RRE
        break;
    else
        i_1=i_2;
    end
end
end
