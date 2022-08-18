function solution=SKM(A,b,xinit,lambda,error_bound,maxIter,beta) 

% xinit is the starting point
% lambda is the projection constant in (0,2]
% error_bound is the upper bound for the relative residual error
% maxIter is the maximum number of iterations 
% beta is the current sample size

%array of iteration indices in which the desired relative residual error is reached
global ITER  


maxIterReached=1;          %variable that indicates if the maximum number of 
	                       %iterations (maxIter) was reached

[m,n] = size(A);			

%column vector with norm of each row
norm_rowsA=sqrt(sum((A).^2,2));

%normalize vector b
b(1:m) = b(1:m)./norm_rowsA;   

%normalize matrix A
A(1:m,:) = A(1:m,:)./repmat(norm_rowsA,1,n); 

initres = A*xinit - b; %initial residual
x = xinit;

%subset of violated constraints in initial residual
initres = initres(initres > 0);    


count=1;
while count <= maxIter
    
    %find the appropriate hyperplane to project
    [l,c] = maxviolatedofsample_constraint(A,b,x,beta);

    % If a selected constraint is violated, project on a direction 
    % orthogonal to violated constraint
    if l > 0
        x_0 = x - lambda*(l/(c*c')) * c';  
        x = x_0;  
    end
    
    
    residual = A*x-b;
    %subset of violated constraints in residual
    residual = residual(residual > 0); 

    
    %if the upper bound for the error is satisfied, stop
    if norm(residual,2)/norm(initres,2) < error_bound 
        ITER=[ITER,count];
        count=maxIter; 	               % stop while loop
        solution=x;
	    maxIterReached=0; 		       % the maximum number of iterations 
				       % was not reached
    end
    count=count+1;
end

if(maxIterReached) %if count reached maxIter
    disp('reached maxIter');
ITER=[ITER,count-1];
solution=x;

end

end
