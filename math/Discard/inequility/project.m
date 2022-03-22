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
    suminc = 0;
     for j=1:m
        Ajx=A(j,:)*x;
        inc = 0;
        if Ajx > b(j)+d
            inc=-(b(j)+d-Ajx);
        end
        if Ajx < b(j)-d
            inc=(b(j)-d-Ajx);
        end
        suminc=suminc+inc*inc;
    end
    arr(1,i)=suminc;
end
end

