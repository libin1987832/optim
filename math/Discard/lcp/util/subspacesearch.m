function x=subspacesearch(x0,M,q)
s=issparse(M);
I=(x0>0);
MII=M(I,I);
qI=q(I);
if s
x=sparse(size(x0,1),1);
[x(I),flag] = pcg(MII,-qI,1e-12,30);
else
x=zeros(size(x0));
x(I)=-MII\qI;
end

