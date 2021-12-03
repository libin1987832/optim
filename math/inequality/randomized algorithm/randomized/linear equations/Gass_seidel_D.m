function [x,normr]= Gass_seidel_D(A, r, maxit,Acol,alpha)

%%
[m, n] = size(A);


%% º∆À„≤–≤Ó
x = zeros(n,1);

for i = 1:maxit
    pickedj=mod(i-1,n)+1;
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;
end
