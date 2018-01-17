function [x,z,Ak,nextA,f,kktz,kkta]=PADSA(H,c,A,n)
x=zeros(n,1);
z=zeros(n,1);
N=1:n;
I=setdiff(N,A);
if ~isempty(A)
    x(I)=-H(I,I)\(c(I)+H(I,A)*x(A));
else
    x(I)=-H(I,I)\(c(I));
end
zt=-1*(H*x+c);
z(A)=zt(A);
Ak=[];
nextA=intmin('uint32');
for i=1:n
    if ismember(i,A) 
        if z(i)>=0;
            Ak=[Ak,i];
            nextA = bitset(nextA,i,1);
        end
    else if x(i)>0
            Ak=[Ak,i];
            nextA = bitset(nextA,i,1);
        end   
    end
end
f=0.5*x'*H*x+c'*x;
kz=z;
kz(kz>0)=0;
kktx=x;
kktx(kktx<0)=0;
kktz=norm(kz);
kkta=norm(kktx);