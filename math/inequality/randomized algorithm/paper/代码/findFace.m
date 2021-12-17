function [index1,index2] = findFace(err,r) 
[m,n]=size(err);
tol = 1e-15;
index1 = -1;
for i=1:n
    rk=err(:,i);
    ssign=sum(~xor(rk>tol, r>tol));
    if ssign == m && index1<0
        index1 = i;
        index2 = i;
    end
    if ssign ~= m && index1>0
        index2 = i;
    end
end