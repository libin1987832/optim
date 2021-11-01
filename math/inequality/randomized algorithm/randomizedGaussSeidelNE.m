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
e = 1;
xA = [x0];
normrow = [];
index = [];
indexA=[0];
%compute norm per row also store the corresponding index
%alpha = 1+ min(m/n,n/m);
alpha = 1;
  for j = 1:n
    normrow = [normrow,norm(A(:,j))];
    index = [index,j];
  end
  
  weight = normrow/sum(normrow);
  Acol=sum(A.*A,1);
if isempty(tol)
  iter = maxit;   
  for i = 1:maxit
    %randsample to generate weighted random number from given vector
    pickedj = randsample(index,1,true,weight);
 %   pickedj=index(mod(i,n)+1);
    indexA = [indexA,pickedj];
    r=b-A*x;
    r(r<0)=0;
     
    col = A(:, pickedj);
    x(pickedj) = x(pickedj) + alpha*( col' * r ) / Acol(pickedj);
      xA =[xA x];
    e = norm(x-exactx);
    error = [error,e];
    %iter = iter+1;
  end
else
    iter = 1;
    e = 1;
    while e >= tol  
    %randsample to generate weighted random number from given vector
    pickedi = randsample(index,1,true,weight);
    
    row = A(pickedi, :);
    x = x + ( b(pickedi) - (row * x) ) / (row * row') * row';
    e = norm(x-exactx);
    error = [error,e];
    iter = iter+1;
    end
end
end