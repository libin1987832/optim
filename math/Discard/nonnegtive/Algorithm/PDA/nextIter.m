function [x,z,Ak,nextA]=nextIter(Ht,b,A,n,x0)
% [x,z]=solver(Ht,b,A,n,x0);
N=1:n;
I=setdiff(N,A);
x=zeros(n,1);


[x,f]=alg2(Ht,b,x0,0.00001);
dig=b-Ht*x;
D=max(diag(dig),0);
D(D>0)=1;
z=Ht'*D*dig;
Ak=[];
nextA=intmin('uint32');
for i=1:n
    if ismember(i,A) 
        if z(i)<=0;
            Ak=[Ak,i];
            nextA = bitset(nextA,i,1);
        end
    else if x(i)<0
            Ak=[Ak,i];
            nextA = bitset(nextA,i,1);
        end   
    end
end
