function [classical_iter_or_error,simple_iter_or_error,rand_iter_or_error] = driver_as_function(iter,tol)
%% set the dimention
r = 10;
m = 100;
n = 2*r+1;

%% generate nodes
%tj are drawing randomly from a uniform distribution in [0,1]
t = rand(m,1);

%then ordering by magnitude (since they are all positive, just sort)
t = sort(t);
% just to assign the size
w = zeros(m,1); 
x = zeros(n,1);
%% generate x 
realx = randn(n,1);
imgx = randn(n,1);
for l = 1:n
   x(l) = realx(l)+1i*imgx(l);
end  
%% generate A and b
A = zeros(m,n);
for j = 1:m
    % dealing with special cases when reach the endpoints(nodes)
    if j == 1
        w(j) = (t(2)-t(end)-1)/2;
    elseif j == m
            w(j) = (t(1)+1-t(m-1))/2;
        else
            w(j) = (t(j+1)-t(j-1))/2;
    end              
    for k = -r:r
       A(j,k+r+1) =sqrt(w(j))*exp(2*pi*1i*k*t(j));
    end
end  
A = 2 * rand(m , n)-1;
b = A*x;
x0 = zeros(n,1);
SP=eye(n^2);
sumA=sum(A.*A,2);
Bn=diag(1./sumA)*A;
for i=1:m
    B=Bn(i,:);
    P=eye(n)-B'*B;
    kp=kron(P,P);
    SP=SP*kp;
end
max(eig(SP))
SP=zeros(n^2,n^2);
for i=1:m
    B=Bn(i,:);
    P=eye(n)-B'*B;
    p=1/m;
    kp=p*kron(P,P);
    SP=SP+kp;
end
max(eig(SP))^m
SP=zeros(n^2,n^2);
normrow = [];
 for i = 1:m
    normrow = [normrow,norm(A(i,:))];
 end  
  weight = normrow/sum(normrow);
for i=1:m
    B=Bn(i,:);
    P=eye(n)-B'*B;
    p= weight(i);
    kp=p*kron(P,P);
    SP=SP+kp;
end
max(eig(SP))^m
[x1,iter1,error1] = kaczmarz(A,b,x0,iter,tol,x);
[x2,iter2,error2] = simplerandkaczmarz(A,b,x0,iter,tol,x);
[x3,iter3,error3] = randkaczmarz(A,b,x0,iter,tol,x);
if isempty(tol)
    classical_iter_or_error = error1;
    simple_iter_or_error = error2;
    rand_iter_or_error = error3;
else
    classical_iter_or_error = iter1;
    simple_iter_or_error = iter2;
    rand_iter_or_error = iter3; 
end