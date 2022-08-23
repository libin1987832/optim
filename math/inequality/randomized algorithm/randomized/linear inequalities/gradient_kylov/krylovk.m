% 用LSQR算法计算子问题 argmin F（u） = || Au + Ax_k-b-z_k || = || Au -y_k ||
function xk=krylovk(A,y,k)
u1=0;

beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;

ro_1=alph1;
thgma_1=beta1;
g1=v1;
[m,n]=size(A);
Av1=zeros(m,1);
Aq2=zeros(n,1);
for i=1:k
    
    for j = 1:m
       Av1(j) = A(j,:)*v1;
    end
%     q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
    q2=Av1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
    ro1=norm([ro_1,beta2]);
    
    c1=ro_1/ro1;s1=beta2/ro1;
    for j = 1:n
       Aq2(j) = A(:,j)'*q2;
    end
   % v2=A'*q2-beta2*v1;alph2=norm(v2);
   v2=Aq2-beta2*v1;alph2=norm(v2);
    v2=v2/alph2;
    
    theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
    
    u2=u1+thgma1*g1./ro1;  g2=v2-theta2*g1./ro1;
    
    u1=u2;
    q1=q2;v1=v2;alph1=alph2;
    
    ro_1=ro_2;
    thgma_1=thgma_2;
    g1=g2;
end
xk=u1;
end