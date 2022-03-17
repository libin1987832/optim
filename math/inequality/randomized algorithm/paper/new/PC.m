function [xk,zk,fk]=PC(x0,z1,A,b,maxit,tol,debug)
[m, n] = size(A);
x = x0;
iter = 0;
if debug
if isempty(exactx)
    r=b-A*x0;
    r(r<0)=0;
    error_k = [norm(A'*r)];
    iter_k =[0];
else
    error_k = [norm(x-exactx)];
    iter_k =[x];
end
end
index_k=[0];
norm_rn = zeros(m,1);
for i = 1:maxit
norm_r=norm_rn;
e1=A'*(A*x1-z1-b);
z=A*x1-b;
z(z<0)=0;
r=-z;
r(r<0)=0;
norm_rn=norm(r);
e2=z1-z;
sum=e1'*e1+e2'*e2;
ae=A*e1-e2;
p=(sum)/(sum+ae'*ae);
xk=x1-p*e1;
zk=z1-p*e2;
   if i>1
        if abs(norm_rn-norm_r)<tol || norm_rn < 1e-6
            break;
        end
   end
end
