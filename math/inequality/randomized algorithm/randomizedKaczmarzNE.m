function [x,iter,error,xA,indexA] = randomizedKaczmarzNE(A, b, x0,maxit,tol,exactx)
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
indexA=[];
normrow = [];
index = [];
%compute norm per row also store the corresponding index

  for i = 1:m
    normrow = [normrow,norm(A(i,:))];
    index = [index,i];
  end
  
  weight = normrow/sum(normrow);
  weightOrig = weight;
  indexOrig = index;
  Arow=sum(A.*A,2);
  update = 0;
if isempty(tol)
  iter = maxit;   
  for i = 1:maxit
    %randsample to generate weighted random number from given vector

    pickedi = randsample(index,1,true,weight);
    indexA = [indexA,pickedi];
    row = A(pickedi, :);
    r=b(pickedi) - (row * x); 
    if r>0
    x = x + ( r ) / (Arow(pickedi)) * row';
    update = update + 1;
    end
    xA =[xA x];
%     if size(index,2) == 2
%          weight = weightOrig;
%          index = indexOrig;
%     end
%     loc=find(index==pickedi);
%     index(loc) = [];
%     weight(loc) = [];
%     weight = weight./sum(weight);

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
    r=b(pickedi) - (row * x); 
    if r>0
        update = update + 1;
        x = x + ( r ) / (Arow(pickedi)) * row';
    elseif size(index,1) == 1
         weight = weightOrig;
         index = indexOrig;
    end
    index(pickedi) = [];
    weight(pickedi) = [];
    weight = weight./sum(weight);
    e = norm(x-exactx);
    error = [error,e];
    iter = iter+1;
    end
end
update
end
