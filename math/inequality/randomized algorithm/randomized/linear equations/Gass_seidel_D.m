function [x,normr]= Gass_seidel_D(A, r, maxit,Acol,alpha)

%%
[m, n] = size(A);


%% º∆À„≤–≤Ó
x = zeros(n,1);

normr=[norm(A'*r)];
for i = 1:maxit
    pickedj=mod(i-1,n)+1;
   % pickedj=pickedj_a(pickedj_i(i));
  %  pickedj=randsample(index,1,true,weight);
    col = A(:, pickedj);
    inc = alpha*( col' * r ) / Acol(pickedj);
  %  inc = alpha*( At_r(pickedj) ) / Acol(pickedj);
    x(pickedj) = x(pickedj) - inc;
    r = r - inc*col;
    normr=[normr,norm(A'*r)];
end
end