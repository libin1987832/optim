function  [xk,flag,relres,iter,resvec,arvec,itersm,tf] =otherAlg(A,b,x0,maxit,type) 
t=clock;

% stop criterion
tol = 1e-13;
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
            for i=1:110
            [xk,zk]=PC(x0,z1,A,b);
            z1=zk;x0=xk;
            end
        case 'H'
            %[xk,rk,countFM,countNW,beginNW,tf,vk,rkArr]=han(x0,A,b,maxit);
             
            I=find(rpk>=tol);
            %��ȡ�Ӿ����ж��Ƿ�����
            AI=A(I,:);
            %    AII=AI'*AI;
             hk=AI\rpk(I);      
            aa=spiecewise(A,b,hk,x0);
            xk=x0+aa*hk;
            rpk=b-A*xk;
        case 'B'
            if iter == 1
                printlevel = -1;
                dirmeth = 1;
                linemeth =2;
                epsline=1e-15;
                neq = 0;
                J=zeros(1,maxit);
                E=zeros(1,1);
            end
            r=rpk;
            I=find(r >= 0);
            fval = 0.5*r(I)'*r(I);
            %  [xk,fval,gradient,err,J,iter] = ineq(-A,x0,-b,tolerances,dirmeth,linemeth,maxit,printlevel,neq);
            if (dirmeth == 0)
                [d,R,error,E] = searchdir(A,I,J(:,iter),r,0,E);
                if (printlevel > 1)
                    disp(['                SD info: search direction found by COF'])
                end  % if (printlevel > 1)
            elseif (dirmeth == 1)
                [d,~,error,E] = searchdir(A,I,J(:,iter),r,1,E);
                if (printlevel > 1)
                    disp(['                SD info: search direction found by SVD'])
                end  % if (printlevel > 1)
            elseif (dirmeth == 2)
                [d,R,error,E] = searchdir(A,I,J(:,iter),r,2,E);
                if (printlevel > 1)
                    disp(['                SD info: search direction found by QRE'])
                end  % if (printlevel > 1)
            elseif (dirmeth == 3)
                if (iter == 1 )
                    [d,R,error,E] = searchdir(A,I,J(:,iter),r,2,E);
                    if (printlevel > 1)
                        disp(['                SD info: search direction found by QRE'])
                    end  % if (printlevel > 1)
                else % dirmeth = 3 and iter > 1 in this case
                    if (m <= n ),
                        [d,R,error,E] = searchdir(A,I,J(:,iter-1),r,2,E);
                        if (printlevel > 1)
                            disp(['                SD info: search direction found by QRE'])
                        end  % if (printlevel > 1)
                    else, % try update/downdate
                        if (printlevel > 1)
                            disp(['                SD info: attempting up/downdate '])
                        end  % if (printlevel > 1)
                        [d,R,error,E] = searchdir(A,I,J(:,iter-1),r,3,E);
                        if (error == 1)
                            [d,R,error,E] = searchdir(A,I,J(:,iter-1),r,2,E);
                            if (printlevel > 1)
                                disp(['                SD info: search direction found by QRE'])
                            end  % if (printlevel > 1)
                        else
                            if (printlevel > 1)
                                disp(['                SD info: search direction found by up/downdate'])
                            end  % if (printlevel > 1)
                        end % if (error == 1)
                    end % if (m <= n )
                end % if (iter == 1 )
            end % if (dirmeth == 1)
            
            if (printlevel > 1)
                disp(['                SD info: time required was ',num2str(timereqd)])
            end  % if (printlevel > 1)
            %
            % D stores all the normalized search directions.  This
            % can be used to determine later if some form of conjugacy
            % exists among the directions.
            %
            %    D = [D,d/norm(d)];
            %
            %
            %---------------------------------------------------------
            % Line search:
            %---------------------------------------------------------
            %
           % [stepsize,fval,linederiv,err] = linesrch(A,x0,d,r,fval,epsline,linemeth,neq,printlevel);
               stepsize=spiecewise(A,b,d,x0);
            %
            % Update x to x + stepsize*d;
            %
            xk = x0 + stepsize*d;
                 rpk=b-A*xk;
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