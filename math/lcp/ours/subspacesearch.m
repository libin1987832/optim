function x=subspacesearch(x0,M,q)
s=issparse(M);
I=(x0>0);
MII=M(I,I);
qI=q(I);
if s
x=sparse(size(x0));
x(I) = pcg(MII,-qI,1e-7,150);
else
x=zeros(size(x0));
x(I)=-MII\qI;
end

