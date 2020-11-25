% Dax hybrid algorithm and r=B*r Nr==N nIter 预测的间隔 eIter 预测未来的多少步
function [xk,rk,countFM,countNW,beginNW,tf,vk,rkArr]=sgradientFM(x0,A,b,nIter,diff,maxIter)
t=clock;
rkArr=zeros(2*maxIter,1);
%compute hybrid uIter
[m,n]=size(A);

%FM need a qr decompose
[Q,R]=qr(A);

rkp=b-A*x0;
r=rkp;
r(r<0)=0;
%condition for terminate
Ar=norm(A'*r);
rn=norm(r);
am=max(max(A));
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
uIndex=0;

countFM=0;
countNW=0;
beginNW=0;
L1=0;
L2=0;
while Ar>delt*rn && rn>delt
%while rn>delt
    %     if uIndex<uIter
    if uIndex<nIter
        countFM=countFM+1;
        uIndex=uIndex+1;
        rkp0=rkp;
        %FM algorithm
        [xk,rkp]=FixedM(x0,Q,R,A,b,rkp);
        rk=rkp;
        rk(rk<0)=0;
    end
    
    if mod(uIndex,nIter)==0 && uIndex >0
        uIndex=0;     
        r0=rkp0;
        r0(r0<0)=0;
        %|| L(xk,zk) ||
        L1=r0'*r0;
        %||L(xk+1,zk+1)||^2
        L2=rk'*rk;
        z0=-rkp0;
        z0(z0<0)=0;
        % -(b-Axk+1)-zk
        LMv=-rkp-z0;
        LM=LMv'*LMv;
        %L(xk,zk)-L(xk+1,zk)
        L1m=L1-LM;
        %L(xk+1,zk)-L(xk+1,zk+1)
        L2m=LM-L2;
        %[L(xk,zk)-L(xk+1,zk)]-[L(xk+1,zk)-L(xk+1,zk+1)]
        LL=L1m-L2m;
        
        % active set r0
        AAz=(rkp0>ee);
        %Au=-(b-Ax_{k+1})+(b-Axk)
        Au=-rkp+rkp0;
        % Au(AAz)
        AuAA=Au(AAz);
        %|| Au(AAz) ||^2
        LLA=AuAA'*AuAA;
        
        % if LLA>rou*LL || LLA<elta
        if abs(LLA-LL)<diff
            %newtonalgorithm
            countNW=countNW+1;
            % record begin newtron type iteratror
            if countNW ==1
                beginNW=countFM;
            end
%             if type ==1
%                 [xk,~]=krylov(A,b,xk,rkp);
%             else
%                 I=find(rkp>=ee);
%                 % 提取子矩阵判断是否正定
%                 AI=A(I,:);
%                 hk=AI\rkp(I);
%                 aa=piecewise(A,b,hk,xk);
%                 xk=xk+aa*hk;
%             end
            [xk,rkp]=sms(A,b,xk,rkp);
            rk=rkp;
%             rkp=b-A*xk;
%             rk=rkp;
            rk(rk<0)=0;
          
%             if rou>0.6
%             nIter=(2^countNW)*nIter;
%             end
            %nIter=(2^countNW)*nIter;
            %             [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
            %             xkArr=[xkArr;[xk',fk,1]];
        end
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
    rkArr(countFM+countNW)=rn;
    % test
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));
%disp(['QQ:',num2str(var)]);
%disp(['hybrid6 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);