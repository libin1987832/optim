function [xk,rk]=solver(Ht,b,A,n,x0)
N=1:n;
I=setdiff(N,A);
x=zeros(n,1);

xk=x0;
while(1)
    dig=Ht*xk-b;
    D=max(diag(dig),0);
    D(D>0)=1;
    H=Ht'*D*Ht;
    c=Ht'*D*b;
    de=det(H(I,I));
    if norm(de)<0.0001
         fprintf('H cannot inverse\n');
        break;
    end
    x(I)=H(I,I)\(c(I));
    xk=x;
    dig1=Ht*xk-b;
    if all(dig==dig1)
        fprintf('solver AT*(Ax-b)+ is same with x0\n');
        break;
    end
    fprintf('solver AT*(Ax-b)+ is different from x0\n');
end

rk=H*xk-c;