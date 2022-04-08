m=5;
n=3;
A = rand(m,n);
b = rand(m,1);
x00 = rand(n,1);
maxiter = 100;
normr=zeros(2,maxiter);
 for i = 1:n
   normrow = [normrow,A(:,i)'*A(:,i)];
   index = [index,i];
end
weight = normrow/sum(normrow);
for k = 1:maxiter
    i = randsample(index,1,true,weight);
    r=b-A*x1;
    r(r<0)=0;
    x1(i)=x1(i)+lambda*A(:,i)'*r/(A(:,i)'*A(:,i));
end
x2=x00;
for k = 1:maxiter
    i = mod(k,n) + 1; 
    r=b-A*x2;
    r(r<0)=0;
    x2(i)=x2(i)+lambda*A(:,i)'*r/(A(:,i)'*A(:,i));
end