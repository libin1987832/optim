% pina hybrid algorithm
function [xk,rk,countFM,error_k,beginNW,tf,vk,rkArr]=han(x0,A,b,maxIter,tol,exactx,debug)
t=clock;

%compute hybrid uIter
[m,n]=size(A);
rkArr=zeros(2*maxIter);
%FM need a qr decompose
r0=b-A*x0;
r0(r0<0)=0;
%condition for terminate
norm_Ar=norm(A'*r0);
norm_r=norm(r0);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;

countFM=0;
countNW=0;
beginNW=0;
if isempty(exactx)
    error_k = [norm_Ar];
    iter_k =[0];
else
    error_k = [norm(x-exactx)];
    iter_k =[x];
end
% AII=AI'*AI;
% bII=AI'*(r(I));
% sparse for svds densy for svd
% [U,S,V]=svds(AII);

if norm_Ar<delt*norm_r || norm_r<delt
    xk=x0;
    rk=r0;
    fk=0.5*(r0'*r0);
    disp('input x is satisfied all constrain!(Ar<delt*rn|| rn<delt)') %ceases execution
end
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de

while 1
    countFM=countFM+1;
    I=find(r0>=tol);
    %提取子矩阵判断是否正定
    AI=A(I,:);
%    AII=AI'*AI;
    hk=AI\r0(I);     
     %   size(AI)
    aa=spiecewise(A,b,hk,x0);
    xk=x0+aa*hk;
    
    rk=b-A*xk;
    rk(rk<0)=0;
    
    
    r0=rk;
    x0=xk;
    Ar=norm(A'*rk);
    norm_rn=norm(rk);
    rkArr(countFM)=norm_rn;
     if Ar<tol || norm_rn < 1e-6

            break;
     end
     norm_r=norm_rn;
     if debug
        if ~isempty(exactx)
            e = norm(x-exactx);
            iter_k =[iter_k x];
        else
            e = Ar;
            iter_k =[iter_k i];
        end
        error_k = [error_k,e];
    end
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk)>0);
%disp(['%hybrid1 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);
%disp(['$',num2str(m),'\times ',num2str(n),'$&FM&(',num2str(countFM),',',num2str(countNW),')&',num2str(tf),'&',num2str(fk),'&',num2str(Ar)]);
%disp(['well1033&Daxs&',num2str(vk),'&',num2str(rn),'&',num2str(Ar),'&(',num2str(countFM),',',num2str(countNW),')&',num2str(beginNW)]);


