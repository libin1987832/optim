function [x,error]= GS(A, b, RRE,x_0,alpha,maxit,b_A,debug)
%%
[~, n] = size(A);
%% º∆À„≤–≤Ó
x = x_0;
Acol=sum(A.*A,1);
error=zeros(1,maxit+1);
r=b-A*x_0;
error(1)=r'*r;
for i = 1:maxit
    pickedj=mod(i-1,n)+1;
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;
     if debug
            error(i+1)=r'*r;
     end
    if  norm(r-b_A)^2<RRE
        break;
    end
end