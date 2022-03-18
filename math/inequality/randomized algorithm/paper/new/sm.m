
% compuation lsqrm: A*x(n)   lsqminnorm: A*x(n) A*xs spiecewise A*x(n)
function [xs,rpk,len,flag]=sm(A,b,n,rpk,x0)
% lsqr tollera
tol = -1e-20;
[m,n] = size(A);
rpk=b-A*x0;
% subspace
for i=1:1
    
    I = find(rpk>=tol);
    AI = A(I,:);
    %       hk = lsqminnorm(AI,rpk(I));
          hk = AI \ rpk(I);
    %         [Q,R]=qr(AI);
    %         qrpk=Q'*rpk(I);
    %         hk=pinv(R)*qrpk;
%     [U,S,V]=svd(AI);
%     minmn = min(size(S));
%     svdj=minmn;
%     for j=1:minmn
%         if S(j,j)<1e-20
%             svdj=j;
%             break;
%         end
%     end
%     B=U(:,1:svdj)'*rpk(I);
%     for j=1:svdj
%         B(j)=B(j)/S(j,j);
%     end
%     hk=V(:,1:svdj)*B;
%     
    
    xss = x0 + hk;
    rpks = b(I) - AI*xss;
    if norm(rpks)<1e-13
        Is = setdiff(1:m, I);
        rpkss = b(Is)-A(Is,:)*xss;
        Ass=A(Is,:)*V(:,svdj+1:n);
        [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(zeros(n-svdj,1),Ass,rpkss,100);
        %     [rk, rkh, dh, gh] = residual(Ass,rpkss,xkh);
        xssv=xss+V(:,svdj+1:n)*xkh;
        [rk2, rkh2, dh2, gh2] = residual(A,b,xssv);
        if dh2<1e-13
            xs=xssv;
            rpk = b - A * xs;
            aa = 10;
            break;
        else
            disp('不在解的范围内')
        end
    end
    
    
    aa = spiecewise(A,b,hk,x0);
    xs = x0 + aa * hk;
    
    x0 = xs;
    rpk = b - A * x0;
    
end
len = aa;
flag = 1 + 6;
end
%rpk = b - A * xs;
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