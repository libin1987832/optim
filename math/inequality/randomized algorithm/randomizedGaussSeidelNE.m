function [x,iter,error,xA,indexA] = randomizedGaussSeidelNE(A, b, x0,maxit,tol,exactx)
% randomized kaczmarz by Algorithm 1
% Ax = b
% A - input matrix
% b - right vector
% x0 - initial x
%
% x - the approximated x
% iter - number of iterations before convergence (= maxit if maxit is given)
% error - the norm of difference in x and exact solution after every
% iteration
m = size(A,1);
n = size(A,2);

x = x0;
%iter = 0;
error = [];
  r=b-A*x;
e = 1;
xA = [x0];
normrow = [];
index = [];
indexA=[0];
%compute norm per row also store the corresponding index
%alpha = 1+ min(m/n,n/m);
alpha = 1;
rp=r;
rp(rp<0)=0;
ATrp = A'*rp;

Acol=sum(A.*A,1);
  for j = 1:n
    normrow = [normrow,sqrt(Acol(j))];
    index = [index,j];
  end
  weight = normrow/sum(normrow);

if isempty(tol)
  iter = maxit;   
  for i = 1:maxit
    %randsample to generate weighted random number from given vector
%     normrow= ATrp.*ATrp;
%     weight = normrow/sum(normrow);
    pickedj = randsample(index,1,true,weight);
 %   pickedj=index(mod(i,n)+1);
    indexA = [indexA,pickedj];

    col = A(:, pickedj);
    inc = alpha*( col' * rp ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    r = r - inc*col;
    rp=r;
    rp(rp<0)=0;
%     ATrp = A'*rp;
    xA =[xA x];
    e = norm(x-exactx);
    error = [error,e];
    %iter = iter+1;
  end
else
        iter =0;
           xA=[0];

      [~, ~, normr, normAr] = residual(A,b,x0);
                        error = [normAr];
    while normAr > tol  && normr > tol
    %randsample to generate weighted random number from given vector
    pickedj = randsample(index,1,true,weight);
 %   pickedj=index(mod(i,n)+1);
    indexA = [indexA,pickedj];

    col = A(:, pickedj);
    inc = alpha*( col' * rp ) / Acol(pickedj);
    x(pickedj) = x(pickedj) + inc;
    r = r - inc*col;
    rp=r;
    rp(rp<0)=0;
%     ATrp = A'*rp;

    iter = iter+1;
    
    if iter > maxit
        return;
    end
    if mod(iter,1000)==0
        [~, ~, normr, normAr] = residual(A,b,x);
        xA =[xA iter];
        error = [error,normAr];
    end
    end
end
end