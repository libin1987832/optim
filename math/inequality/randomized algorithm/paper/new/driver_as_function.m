function [classical_iter_or_error,simple_iter_or_error,rand_iter_or_error] = driver_as_function(A,b,x0,r_exact,maxit,tol)
%% set the dimention
% r = 1;
% m = 100;
% n = 2*r+1;

% % generate nodes
% tj are drawing randomly from a uniform distribution in [0,1]
% t = rand(m,1);
% 
% then ordering by magnitude (since they are all positive, just sort)
% t = sort(t);
% just to assign the size
% w = zeros(m,1); 
% x = zeros(n,1);
% % generate x 
% realx = randn(n,1);
% imgx = randn(n,1);
% for l = 1:n
%    x(l) = realx(l)+1i*imgx(l);
% end  
% % generate A and b
% A = zeros(m,n);
% for j = 1:m
%     dealing with special cases when reach the endpoints(nodes)
%     if j == 1
%         w(j) = (t(2)-t(end)-1)/2;
%     elseif j == m
%             w(j) = (t(1)+1-t(m-1))/2;
%         else
%             w(j) = (t(j+1)-t(j-1))/2;
%     end              
%     for k = -r:r
%        A(j,k+r+1) =sqrt(w(j))*exp(2*pi*1i*k*t(j));
%     end
% end  
% A = 2 * rand(m , n)-1;
% b = A*x;
% x0 = zeros(n,1);
[m,n]=size(A);

[x1,iter1,error1,~,~]= GuassSeidelNE(A, b, x0,2 ,maxit,-1,r_exact,1);
[x2,iter2,error2,~,~]=  simpleGuassSeidelNE(A, b, x0,2 ,maxit,-1,r_exact,1);
[x3,iter3,error3,~,~] = randGuassSeidelNE(A, b, x0,2 ,maxit,-1,r_exact,1);
if isempty(tol)
    classical_iter_or_error = error1;
    simple_iter_or_error = error2;
    rand_iter_or_error = error3;
else
    classical_iter_or_error = iter1;
    simple_iter_or_error = iter2;
    rand_iter_or_error = iter3; 
end