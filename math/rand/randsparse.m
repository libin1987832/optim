function [X,IT,errornorm,Ynorm]=randsparse(A,b,y,u,K)
[m,n]=size(A);

if nargin < 4   
    error('Too few input arguments')  
end
if nargin < 5    
    K=10000;
end
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
end

number=zeros(1,m);
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
x=zeros(n,1);
z=zeros(n,1);
k=0;
ynorm1=norm(y);
ynorm=norm(y-x)/norm(y);

while ynorm>=1e-3&&k<K
%while norm(A*x-b,2)>=1e-3&&k<=K
    k=k+1;
    X(:,k)=x;
    IT(k)=k;                                  
    errornorm(k)=norm(A*x-b);
    Ynorm(k)=ynorm;
    i=randsrc(1,1,[number;prob]);
    norm2=prob(i)/F;
    z=z+(b(i)-A(i,:)*x)/norm2*A(i,:)'; 
    for j=1:n
        if z(j)>u
            x(j)=z(j)-u;
        elseif z(j)<-u
            x(j)=z(j)+u;
        else
            x(j)=0;
        end
    end
    ynorm=norm(y-x)/ynorm1;
end
IT=IT(:,1:k);
X = X(:,1:k);
errornorm=errornorm(1,1:k);
Ynorm=Ynorm(1,1:k);

