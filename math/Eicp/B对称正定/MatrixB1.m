function B1=MatrixB1(n)
B1=diag(ones(n,1)*10)+diag(ones(n-1,1)*(-4),-1)+diag(ones(n-1,1)*(-4),1)+diag(ones(n-2,1),2)+diag(ones(n-2,1),-2);
end