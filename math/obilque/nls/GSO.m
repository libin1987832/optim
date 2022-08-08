function [x,error]=GSO(A,b,RRE,x0,tol,iter,b_A,debug)
[~,n]=size(A);
%b_A=b-(A'\(A'*b));
RRE=RRE*norm(b,2)^(2);
error=zeros(1,iter+1);
%N=zeros(n,1);%���ڴ洢�����з���
% for i=1:1:n
%     N(i)=norm(A(:,i),2)^(2);
% end
N=sum(A.*A,1);
%��һ�ε���
r=b-A*x0;
rp=r;
rp(rp<0)=0;
if debug
    error(1)=rp'*rp;
end
x=x0;
alpha=A(:,1)'*rp/N(1);
x(1)=x(1)+alpha;
r=r-alpha*A(:,1);
rp=r;
rp(rp<0)=0;
i_1=1;
for k=1:1:iter
    i_2=mod(k,n)+1;
    G=A(:,i_1)'*A(:,i_2);
    g=N(i_2)-G^2/N(i_1);
    if g>tol
        alpha=A(:,i_2)'*rp/g;
        beta=-G*alpha/N(i_1);
        ab=[alpha;beta];
        x([i_2,i_1],1)=x([i_2,i_1],1)+ab;
        r=r-[A(:,i_2),A(:,i_1)]*ab;
        rp=r;
        rp(rp<0)=0;
        if debug
            error(k+1)=rp'*rp;
        end
    end
    if  norm(r-b_A)^2<RRE
        break;
    else
        i_1=i_2;
    end
end
end
