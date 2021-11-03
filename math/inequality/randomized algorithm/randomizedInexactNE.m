function [x,iter,error,xA,indexA] = randomizedInexactNE(A, b, x0,maxit,tol,exactx)
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
t0=1;
r=b-A*x0;
z=-r;
z(z<0)=0;
z0=z;
alpha = 1+ min(m/n,n/m);
  Acol=sum(A.*A,1);
  for j = 1:n
      normrow = [normrow,sqrt(Acol(j))];
     index = [index,j];
  end
  weight = normrow/sum(normrow);
if isempty(tol)
  iter = 0;   
  for i = 1:maxit
    %randsample to generate weighted random number from given vector
     pickedj = randsample(index,1,true,weight);
  %  pickedj=index(mod(i,n)+1);
     indexA = [indexA,pickedj];
     col = A(:, pickedj);
     rz=(r+z);

     inc = alpha*( col' * rz) / Acol(pickedj);
     x(pickedj) = x(pickedj) + inc;
     r = r - inc*col;
     
     z=-r;
     z(z<0)=0;
     t1 = 0.5 + 0.5 * sqrt(1+4*t0^2);
     z=z0+t0/t1*(z-z0);
     z0=z;
     t0=t1;
     
    xA =[xA x];
    e = norm(x-exactx);
    error = [error,e];
    iter =iter+1;
    %iter = iter+1;
  end
else
    [~, ~, normr, normAr] = residual(A,b,x0);
    iter =0;
       xA=[0];
       error = [normAr];
    while normAr > tol  && normr > tol
    %randsample to generate weighted random number from given vector
      pickedj = randsample(index,1,true,weight);
  %  pickedj=index(mod(i,n)+1);
     indexA = [indexA,pickedj];
     col = A(:, pickedj);
     rz=(r+z);

     inc = alpha*( col' * rz) / Acol(pickedj);
     x(pickedj) = x(pickedj) + inc;
     r = r - inc*col;
     
     z=-r;
     z(z<0)=0;
     t1 = 0.5 + 0.5 * sqrt(1+4*t0^2);
     z=z0+t0/t1*(z-z0);
     z0=z;
     t0=t1;
     
   
  
    iter = iter+1;
    if iter > maxit
        return;
    end
    if mod(iter,1000)==0
        [~, ~, normr, normAr] = residual(A,b,x);
         xA =[xA iter];
        error = [error normAr];
    end
end
end