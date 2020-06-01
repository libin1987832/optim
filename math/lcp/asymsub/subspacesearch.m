function x=subspacesearch(x0,M,q)
I=(x0>0);
MII=M(I,I);
qI=q(I);
x=zeros(size(x0));
x(I)=MII\qI;

