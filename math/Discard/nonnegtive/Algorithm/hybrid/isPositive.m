function isTF=isPositive(a)
at = a;
[mt,nt] = size(a,1)
ai = zeros([mt,nt])
ma=mt
while ma>1
    ma = ma-1;
    if at(1,1)>0
        m=1/a(1,1);
        b=a(2:n,1);
        c=a(1,2:n);
        d=a(2:n,2:n);
        a=d-c*m*b;
        ai[]
    elseif at(1,1) == 0
        for temp in a
            
    end
    
end
if a(1,1)>0 && size(a,1)==1
    isTF = 1;
else
    isTF = 0;
end
end

tf = issymmetric(A)
d = eig(A)
isposdef = all(d > tol)
issemidef = all(d > -tol)