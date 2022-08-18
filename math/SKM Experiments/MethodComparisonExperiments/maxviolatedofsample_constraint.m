function [l,a] = maxviolatedofsample_constraint(A,b,x,beta) 

%l : amount of violation
%a : row of A (hyperplane normal) that defines the most violated constraint inside 
    %the sample of size beta

m = size(A,1);  %number of inequalities

perm=randsample(m,beta);

[l,ind]=max(A(perm(1:size(perm,1)),:)*x-b(perm(1:size(perm,1))));
a=A(perm(ind),:);

l = max(l,0);

end
