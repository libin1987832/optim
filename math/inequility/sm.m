function [xs,rpk,len,flag]=sm(A,b,n,rpk,x0)
AA=(rpk>0);
AI=A(AA,:);
bI=rpk(AA);
[u,flag,relres,iter,resvec,lsvec,out] = lsqrm(AI,bI,1e-5,5,[],[],zeros(n,1),A,b,x0,AA);
if ~out 
    xs = x0 + u;
    len = iter;
else
 %   u = AI\bI;
    aa = spiecewise(A,b,u,x0);
    xs = x0 + aa * u;
    len = aa;
    flag = flag + 6;
end
rpk = b - A * xs;
%semilogy(1:iter,resvec(1:iter),'b.');
% if flag == 0 && ~out
%     flag = 0;
% else
%     flag = flag+10;
% end


% normr0 = norm(bI);
% rk=rpk;
% rk(rk<0)=0;
% normrk = norm(rk);
% 
% fprintf('sm:subspace:%d iter:%d len:%g r0:%g rk:%g u:%g \n',sum(AA),iter,len,normr0,normrk,norm(u));