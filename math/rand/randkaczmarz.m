function [X,IT,errornorm,Ynorm]=randkaczmarz(A,b,y,K,x0)
[m,n]=size(A);

if nargin < 4   
    K=10000;
    x0 = zeros(n,1);
end
if nargin < 5 || isempty(x0)   
    x0 = zeros(n,1);
end

if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
elseif size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem')
end

number=1:m;
prob=zeros(1,m);
F=1/(norm(A,'fro'))^2;
for i=1:m
    number(i)=i;
    prob(i)=norm(A(i,:))^2*F;
end


IT=zeros(1,K);
X=zeros(n,K);
errornorm=zeros(1,K);
Ynorm=zeros(1,K);
x=x0; 
k=0;
ynorm1=norm(y);
ynorm=norm(y-x)/norm(y);

while ynorm>=1e-3&&k<K
%while norm(A*x-b,2)>=1e-3&&k<=K
    k=k+1;
    X(:,k)=x;
    IT(k)=k;
    errornorm(k)= norm(A*x-b);
    Ynorm(k)=ynorm;
    i=randsrc(1,1,[number;prob]);
    x=x+(b(i)-A(i,:)*x)/(prob(i)/F)*A(i,:)'; 
    ynorm=norm(y-x)/ynorm1;
end
IT=IT(:,1:k);
X = X(:,1:k);
errornorm=errornorm(1,1:k);
Ynorm=Ynorm(1,1:k);

