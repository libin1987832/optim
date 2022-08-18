function [solution,iter] = randkaczmarz(A,b,xinit,error_bound,maxIter) 
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

x = xinit;


normrow = [];
index = [];
%compute norm per row also store the corresponding index

  for i = 1:m
    normrow = [normrow,norm(A(i,:))];
    index = [index,i];
  end
  
  weight = normrow/sum(normrow);
count=1;
while count <= maxIter
    %randsample to generate weighted random number from given vector
    pickedi = randsample(index,1,true,weight);
    row = A(pickedi, :);
    x = x + ( b(pickedi) - (row * x) ) / (row * row') * row';
    residual = A*x-b;
    if norm(residual,Inf)< error_bound 
        iter = count;
        count=maxIter; 	               % stop while loop
        solution=x;
	    maxIterReached=0; 		       % the maximum number of iterations 
				       % was not reached
    end
    count=count+1;
end

if(maxIterReached) %if count reached maxIter
    disp('reached maxIter');
    iter = count -1;
    solution=x;
end
end