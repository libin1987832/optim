function  [xk,flag,relres,iter,resvec,arvec,itersm,tf] =otherAlg(A,b,x0,maxit,type) 
t=clock;

% stop criterion
tol = 1e-12;
[m,n] = size(A);

[rpk, r0, normr, normAr] = residual(A,b,x0);

iter = 0;
% residual vector
resvec = zeros(1,maxit + 1);
% the normal gradient 
arvec = zeros(1,maxit + 1);
% subspace minization
itersm = zeros(1,maxit + 1);
resvec(1) = normr;
arvec(1) = normAr;
indexsm = 0;
% flag 0-4 return lsqr flag
flag = 5;

%while normAr > tol * normA * normr && normr > tol
while normAr > tol  && normr > tol
    iter = iter + 1;
    switch upper(type)
        case 'F'
            if iter == 1
                [Q,R]=qr(A);
            end
            [xk,rpk]=FMQR(x0,Q,R,A,b,rpk);
        case 'I'
              r0=rpk;
              r0(r0<0)=0;
              uk=krylovk(A,r0,3);
              xk=x0+uk;
              rpk=b-A*xk;
        case 'P'
            if iter == 1
                z1 = -rpk;
                z1(z1<0)=0;
            end
            [xk,zk]=PC(x0,z1,A,b);
            z1=zk;
        case 'H'
            %[xk,rk,countFM,countNW,beginNW,tf,vk,rkArr]=han(x0,A,b,maxit);
             
            I=find(rpk>=tol);
            %提取子矩阵判断是否正定
            AI=A(I,:);
            %    AII=AI'*AI;
             hk=AI\rpk(I);      
            aa=spiecewise(A,b,hk,x0);
            xk=x0+aa*hk;
            rpk=b-A*xk;
        case 'B'
            if iter == 1
                      epsout = 1.0E-12;
                        epsline = 1.0E-15;
                         tolerances = [epsout;epsline];
                 printlevel = 2;
                 dirmeth = 2;
                 linemeth =2;
                       neq = 0;
            end
               [xk,fval,gradient,err,J,iters] = ineq(A,x0,b,tolerances,dirmeth,linemeth,maxit,printlevel,neq);
               break;
    otherwise
        error(message('MATLAB:InvalidOp'))
    end
    x0=xk;
    [rpk, r, normr, normAr] = residual(A,b,x0,rpk);
    % record the value of objection function
    resvec(iter + 1) = normr;
    % record the value of the gradient function
    arvec(iter + 1) = normAr;
    if iter > maxit-1 || flag == 0
        break;
    end
    if flag < 5 || flag > 6 
%         if flag < 12
%             disp(['flag:',num2str(flag)]);
%         end
         flag = 5;
    end
end
xk = x0;
relres = normr;
resvec = resvec(1:iter + 1);
itersm = itersm(1:iter + 1);
arvec = arvec(1:iter + 1);
tf = etime(clock,t);