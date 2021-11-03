function [x,iter,error,xA,indexA] = IFM(A, b, x0,maxit,tol,exactx)

m = size(A,1);
n = size(A,2);

x = x0;
%iter = 0;
error = [];
xA = [x0];
indexA=[];
%compute norm per row also store the corresponding index
rk=b-A*x;

k=3;
iter = 0;
if isempty(tol)
%     iter = maxit;
    for i = 1:maxit
        %randsample to generate weighted random number from given vector

        rk(rk<0)=0;
        uk=krylovk(A,rk,k);
        x=x+uk;
        rk = b - A * x;
        e = norm(x-exactx);
        error = [error,e];
        xA =[xA x];

        iter = iter+1;
    end
else
    [~, ~, normr, normAr] = residual(A,b,x0);
        iter =0;
        xA=[0];
        error = [normAr];
    while normAr > tol  && normr > tol
        rk(rk<0)=0;
        uk=krylovk(A,rk,k);
        x=x+uk;
        rk = b - A * x;
        iter = iter+1;
    if iter > maxit
        return;
    end
     [~, ~, normr, normAr] = residual(A,b,x);
        xA =[xA iter];
        error = [error,normAr];
    end
end

function xk=krylovk(A,y,k)
u1=0;

beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;

ro_1=alph1;
thgma_1=beta1;
g1=v1;

for i=1:k
    q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
    
    ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
    
    v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
    
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

end