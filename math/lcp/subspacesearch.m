function d=subspacesearch(x0,M,q)
I=(x0>0);
MII=M(I,I);
qI=q(I);
d=zeros(size(x0));
d(I)=MII\qI-x0(I);

