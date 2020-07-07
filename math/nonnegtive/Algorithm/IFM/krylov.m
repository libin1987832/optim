function [xk,zk]=krylov(AALL,b,x0)
[m,n]=size(AALL);
r0=AALL*x0-b;
z0=r0;
z0(z0<0)=0;
ee=1e-15;% computer floating point arithmetic
AA=(z0<ee);
%FF=setdiff(1:m,AA)';

A=AALL(AA,:);
FF=setdiff(1:m,AA)';
AFF=AALL(FF,:);
%|| A(AA,:)u-(b(AA)-A(AA,:)x) ||
y=b(AA)-A*x0;
% \fi
rmrk=AALL*x0-b-z0;
%watch1=0.5*rmrk'*rmrk;

u1=0;  
beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;

ro_1=alph1;
thgma_1=beta1;
g1=v1;
xk=x0;
fmk=0;
for i=1:n
    q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
    
    ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
    
    v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
    
    theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
    
    d=thgma1*g1./ro1;
    %g2=v2-theta2*g1./ro1;
    u2=u1+d;  g2=v2-theta2*g1./ro1;
    xk=xk+d;
   % xk=x0+u2;
    zkp=fmk;
    fmk=AFF*xk-b(FF);
    %zk=fmk;
    %zk(zk<0)=0;
    
    
    %rmrk=AALL*xk-b-zk;
    %watch2=0.5*rmrk'*rmrk;

    %AAk=(zk<ee);
    %empty=isempty(setdiff(AA,AAk));
    empty=isempty(find(fmk<0,1));
    if ~empty
        x1=xk-d;
        %zk=zkp;
        %zk(zk<0)=0;
        
        Ad=AFF*d;
        z0Ad=zkp./Ad;
        z0Ad(z0Ad<=0)=inf;
        a2=min(z0Ad);
        a2=min(a2,1);
        xk=x1+a2*d;
        zk=z0;
        zk(FF)=zk(FF)-a2*Ad;
        
        %rmrk=AALL*xk-b-zk;
        %watch3=0.5*rmrk'*rmrk;

        break;
    end
    
    u1=u2;
    q1=q2;v1=v2;alph1=alph2;
    
    ro_1=ro_2;
    thgma_1=thgma_2;
    g1=g2;
end


% function [xk,fk]=kyrlov(A,y,k)
% u1=0;  
% beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;
% 
% ro_1=alph1;
% thgma_1=beta1;
% g1=v1;
% 
% 
% for i=1:k
%     q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
%     
%     ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
%     
%     v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
%     
%     theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
%     
%     u2=u1+thgma1*g1./ro1;  g2=v2-theta2*g1./ro1;
%     
%     u1=u2;
%     q1=q2;v1=v2;alph1=alph2;
%     
%     ro_1=ro_2;
%     thgma_1=thgma_2;
%     g1=g2;
% end
%  xk=u1;
%  fk=0.5*(A*xk-y)'*(A*xk-y);
% 
    