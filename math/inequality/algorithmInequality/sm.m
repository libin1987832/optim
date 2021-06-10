
% compuation lsqrm: A*x(n)   lsqminnorm: A*x(n) A*xs spiecewise A*x(n)
function [xs,rpk,len,flag]=sm(A,b,n,rpk,x0)
% lsqr tollera
tol = 1e-20;
[m,n] = size(A);
    rpk=b-A*x0;
% subspace
for i=1:100

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
    aas = spiecewise(A,b,hk,x0);
    ap=A*hk;
    ai=rpk./ap;
    as=sort(ai(ai>0));
    tas=[0;as];
    aa = tas(2)+1e-13;
    aaa = tas(2);
    xs = x0 + aa * hk;
    xsaa=x0 + aaa * hk;
    rpksssa = b - A * xsaa;
    rpksssa(rpksssa<0)=0;
    rpksss = b - A * xs;
    rpkkkk=rpk;
    rpkkkk(rpkkkk<0)=0;
    rpksss(rpksss<0)=0;
    fprintf('rpk=%g,rpks=%g,rpksa=%g,len=%g,(%d,%d)\n',norm(rpkkkk),norm(rpksss),norm( rpksssa),aa , sum(rpk>=tol),sum(rpksss>=tol));
    
    xss = x0 + aas*hk;
    rpks = b(I) - AI*xss;
    Is = setdiff(1:m, I); 
    rpkss = b(Is)-A(Is,:)*xss;
    Isc=rpkss>0;
    rpkss(rpkss<0)=0;
    fprintf('rpk=%g,rpks=%g,len=%g,%d\n',norm(rpkss),norm(rpks),aas,sum(Isc));
       
    x0 = xs;
    rpk = b - A * x0;
    if norm(rpksss)<1e-14
        break;
    end
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