function [xk,rpk]=sms(A,b,x0,rpk0)
% [m,n]=size(A);
% AA=(rpk0>0);
% AI=A(AA,:);
% bI=rpk0(AA);
%[u,flag,relres,iter,resvec,lsvec] = slsqrM(AI,bI,1e-12,n,[],[],zeros(n,1),A,b,x0,AA);
%[u,flag,relres,iter,resvec,lsvec] = lsqr(AI,b,1e-12,n);
                
% if flag == 5 || flag==1
%     xk = x0 + u;
% else
%     aa = spiecewise(A,b,u,x0);
%     xk = x0 + aa * u;
% end
[xk,~]=skrylov(A,b,x0,rpk0);

rpk = ( b - A * xk);

% rk0=rpk0;
% rk0(rk0<0)=0;
% rk=rpk;
% rk(rk<0)=0;
% rk;