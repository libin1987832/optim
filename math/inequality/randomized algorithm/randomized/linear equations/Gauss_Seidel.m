A=[3,4;5,6];
b=[2;3];
[m,n]=size(A);
B=A'*A;
Ab=A'*b;
Acol=diag(B);

alpha =1;
maxit = 2;
%

%%
[m, n] = size(A);


%% º∆À„≤–≤Ó
r = b;
x = zeros(n,1);
At_r = A'*r;

for i = 1:maxit
    pickedj=mod(i-1,n)+1;
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    col = A(:, pickedj);
   % inc = alpha*( col' * r ) / Acol(pickedj);
    inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;

    r = r - inc*col;
        At_r = At_r - inc*B(:,pickedj);
    col'*r
end
