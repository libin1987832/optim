% 
% A=[-1,0;1,0;0,1];
% b=[1;1;2];
% x0=[5;0];
% [xk,fk,f0]=ssqr2(x0,A,b)
% [xk,fk,f0]=ssqr2(xk,A,b)
%A=[2,1;3,1];
A=[0.5,0.3;1,0.7;1,3;0.3,0.6];
b=[1;-1;1;-1];
x0=[0;0];
addpath('FM');
[m,n]=size(A);
[xk1,fk1,xkArr1]=hybrid1(x0,A,b);
[xk2,fk2,xkArr2]=hybrid2(x0,A,b);

% for ii=1:size(xkArr1,1)
%     xkA1=xkArr1(ii,:);
%     xkA2=xkArr2(ii,:);
%     diffx=max(xkA1(1:n)-xkA2(1:n));
%     difff=xkA1(n+1)-xkA2(n+1);
%     diffm=xkA1(n+2)-xkA2(n+2);
% end
    