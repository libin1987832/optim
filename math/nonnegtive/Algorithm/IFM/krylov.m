function [xk,zk,rkp]=krylov(AALL,b,x0,rkp0)
[m,n]=size(AALL);
z0=-rkp0;
z0(z0<0)=0;
ee=1e-15;% computer floating point arithmetic
AA=(z0<ee);
%FF=setdiff(1:m,AA)';
zk=z0;
A=AALL(AA,:);
FF=~AA;
AFF=AALL(FF,:);
%|| A(AA,:)u-(b(AA)-A(AA,:)x) ||
y=rkp0(AA);
% \fi
%rmrk=AALL*x0-b-z0;
%watch1=0.5*rmrk'*rmrk;

u1=0;  
beta1=norm(y);q1=y/beta1;v1=A'*q1;alph1=norm(v1);v1=v1/alph1;

ro_1=alph1;
thgma_1=beta1;
g1=v1;
xk=x0;
fmk=0;
rkp=rkp0;
for i=1:n
    q2=A*v1-alph1*q1;beta2=norm(q2);q2=q2/beta2;
    
    ro1=norm([ro_1,beta2]);c1=ro_1/ro1;s1=beta2/ro1;
    
    v2=A'*q2-beta2*v1;alph2=norm(v2);v2=v2/alph2;
    
    theta2=s1*alph2;ro_2=c1*alph2;thgma1=c1*thgma_1;thgma_2=-1*s1*thgma_1;
    
    d=thgma1*g1./ro1;
    %g2=v2-theta2*g1./ro1;
    u2=u1+d;  g2=v2-theta2*g1./ro1;
    
    
  %  xk=xk+u2;
  %check I(xk)==I(xk+1)
    xk=x0+u2;
    rkp0=rkp;
    rkp=b-AALL*xk;
    AAk=(rkp>-ee);
    empty=sum(xor(AA,AAk));
    if empty
        if u1==0
            for j=0:5
            I=find(rkp0>=ee);
            AI=AALL(I,:);
            hk=AI\rkp0(I);
            aa=piecewise(AALL,b,hk,x0);
            xk=x0+aa*hk;
          
%             rk0=rkp0;
%             rk0(rk0<0)=0;
%             gradient=AALL'*rk0;
%             aa=piecewise(AALL,b,gradient,x0);
%             xk=x0+aa*gradient;

            x0=xk;
            rkp0=b-AALL*x0;
           end
%             rkp=b-AALL*xk;
%             rk=rkp;
%             rk(rk<0)=0;
%             rr=norm(rk);
%      [xk,rk,countFM,countNW,beginNW,tf,vk,xkArr]=han(x0,AALL,b,0);
% rkp=b-AALL*xk;
        else
            xk=x0+u1;
        end
        rkp=rkp0;
        break;
    end
%     zkp=fmk;
%     fmk=AFF*xk-b(FF);
%     zk(FF)=fmk;
    %zk(zk<0)=0;
    
    
    %rmrk=AALL*xk-b-zk;
    %watch2=0.5*rmrk'*rmrk;

    %AAk=(zk<ee);
    %empty=isempty(setdiff(AA,AAk));
%     empty=isempty(find(fmk<0,1));
%     if empty
%         x1=xk-d;
%         %zk=zkp;
%         %zk(zk<0)=0;
%         
%         Ad=AFF*d;
%         z0Ad=zkp./Ad;
%         z0Ad(z0Ad<=0)=inf;
%         a2=min(z0Ad);
%         a2=min(a2,1);
%         xk=x1+a2*d;
%         zk=z0;
%         zk(FF)=zk(FF)-a2*Ad;
%         
%         %rmrk=AALL*xk-b-zk;
%         %watch3=0.5*rmrk'*rmrk;
% 
%         break;
%     end
    % update paramters
    u1=u2;
    q1=q2;v1=v2;alph1=alph2;
    
    ro_1=ro_2;
    thgma_1=thgma_2;
    g1=g2;
    
    zk=0;
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
    