A=rand(3,3);
b=rand(3,1);
x00=rand(3,1);
x0=x1;
for i=1:3
    x1=x0+A(:,i)'*(b-A*x0)/(A(:,i)'*A(:,i));
    x0=x1;
end
ATA=A'*A;
x01=x