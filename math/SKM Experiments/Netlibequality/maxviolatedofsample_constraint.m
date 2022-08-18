function [l,a] = maxviolatedofsample_constraint(A,b,x) 

%l : amount of violation
%a : row of A (hyperplane normal) that defines the most violated constraint inside 
    %the sample of size beta

m = size(A,1);  %number of inequalities

[l,ind]=max(A*x-b);
a=A(ind,:);

l = max(l,0);

end
