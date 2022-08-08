function [x,arr]=project(A,b,x0,iter,d,l,u,debug)
[m,n]=size(A);
x=x0;
ATA=sum(A.*A,2);
arr=zeros(1,iter);
for i=1:iter
    for j=1:m
        Ajx=A(j,:)*x;
        inc = 0;
        if Ajx > b(j)+d
            inc=(b(j)+d-Ajx)/ATA(j);
        end
        if Ajx < b(j)-d
            inc=(b(j)-d-Ajx)/ATA(j);
        end
        x=x+inc*A(j,:)';
    end
x=max(l,x);
x=min(u,x);
if debug
    Ax=A*x;
    r1=b-Ax-d;
    r1(r1<0)=0;
    r2=-b-d+Ax;
    r2(r2<0)=0;
    arr(1,i)=norm([r1; r2]);
end
end

