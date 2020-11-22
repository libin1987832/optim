function xk = gauss_siedel(xk,A,b,niter)
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;
e= max(eig(-inv(D+L)*(U)));
if abs(e) >= 1
    disp ('Since the modulus of the largest Eigen value of iterative matrix is not less than 1') 
    return
end
for i = 1:niter
    xk = -inv(D+L)*(U)*xk + inv(D+L)*b;% Gauss-Seidel formul   
end
