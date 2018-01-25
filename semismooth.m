function [y,r,h]=semismooth(H,c,y0,r0,n)
A1=[H,diag(ones(n,1))];
A2=zeros(n,2*n);
diff=y0+r0;
delt=0;
for i=1:n
    if diff(i)>0
        A2(i,i)=-1;
    end
    if diff(i)<0
        A2(i,i+n)=1;
    end
    if diff(i)==0
        A2(i,i)=-1*delt;
        A2(i,i+n)=1-delt;
    end
end
F1=H*y0+c+r0;
F2=r0-max(0,diff);
h=-1*[A1;A2]\[F1;F2];
y=y0+h(1:n);
r=r0+h(n+1:2*n);
[A1;A2];

% E=diag(n,n);
% Z=zeros(n,n);
% HH=[H(I,I),H(I,A),E(I,I),Z(I,A);
%     H(A,I),H(A,A),Z(A,I),E(A,A);
%     Z(I,I),Z(I,A),E(I,I),Z(I,A);
%     Z(A,I),-1*E(A,A),Z(A,I),z(A,A)];
% cc1=A*y0+r0-c;
% cc2=r0;
% cc=-1*[H*y0+]
% AII=H(I,I);
% AIA=H(I,A);
% EII=E(I,I);
% AAI=H(A,I);
% AAA=H();