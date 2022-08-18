function [solution,iter]=skmeq(A,b,xinit,error_bound,maxIter) 

% xinit is the starting point
% lambda is the projection constant in (0,2]
% error_bound is the upper bound for the relative residual error
% maxIter is the maximum number of iterations 
% beta is the current sample size

%array of iteration indices in which the desired relative residual error is reached 


maxIterReached=1;          %variable that indicates if the maximum number of 
	                       %iterations (maxIter) was reached

[m,n] = size(A);			

%column vector with norm of each row
norm_rowsA=sqrt(sum((A).^2,2));

%normalize vector b
b(1:m) = b(1:m)./norm_rowsA;   

%normalize matrix A
A(1:m,:) = A(1:m,:)./repmat(norm_rowsA,1,n); 

x = xinit;

%subset of violated constraints in initial residual  


count=1;
while count <= maxIter
    
    %find the appropriate hyperplane to project
    [l,c] = maxviolatedofsample_constraint(A,b,x);

    x_0 = x - l * c';
    x = x_0;  
    
    
    residual = A*x-b;
 

    %if the upper bound for the error is satisfied, stop
    if norm(residual,Inf)< error_bound 
        iter = count;
        count=maxIter; 	               % stop while loop
        solution=x;
	    maxIterReached=0; 		       % the maximum number of iterations 
				       % was not reached
    end
    count=count+1;
end

if(maxIterReached) %if count reached maxIter
    disp('reached maxIter');
    iter = count -1;
    solution=x;
end

end
