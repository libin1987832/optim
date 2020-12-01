function [X,IT,errornorm,Ynorm]=randkaczmarz(A,b,y,maxit,x0)
[m,n]=size(A);

if nargin < 4   
    maxit=1e4;
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
As = A.*A;
F=sum(sum(As));
prob = sum(As,2)/F;



IT=zeros(1,maxit);
X=zeros(n,maxit);
errornorm=zeros(1,maxit);
Ynorm=zeros(1,maxit);
x=x0; 
k=0;
ynorm1=norm(y);
ynorm=norm(y-x)/norm(y);

while ynorm>=1e-3 && k < maxit
%while norm(A*x-b,2)>=1e-3&&k<=K
    k=k+1;
    X(:,k)=x;
    IT(k)=k;
    errornorm(k)= norm(A*x-b);
    Ynorm(k)=ynorm;
    i=randsrc(1,1,[number;prob]);
    t = (b(i)-A(i,:)*x)/sum(As(i,:));
    x=x+ t * A(i,:)'; 
    ynorm=norm(y-x)/ynorm1;
end
IT = IT(:,1:maxit);
X = X(:,1:maxit);
errornorm=errornorm(1,1:maxit);
Ynorm=Ynorm(1,1:maxit);

