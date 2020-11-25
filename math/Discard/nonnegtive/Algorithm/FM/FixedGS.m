% x0 initial value  b
% xk out value 
%r0=(b-Ax0) rk=(b-Axk)  
function [xk,rk]=FixedM_GS(x0,A,b,r,niter)
[m,n]=size(A);
z=-r;
z(z<0)=0;
bz=b+z;
xk=x0;
% compute min increase
for j=1:niter
    for i=1:n
        zk=A*xk-bz;
        c=A(:,i);
        xk(i)=xk(i)-(c'*zk)/(c'*c);
    end
end
rk=b-A*xk;
