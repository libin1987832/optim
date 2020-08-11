% Dax hybrid algorithm and r=B*r Nr==N nIter 预测的间隔 eIter 预测未来的多少步
function [xk,rk,countFM,countNW,beginNW,tf,vk]=gradientFM_i(x0,A,b,nIter,diff,maxIter,xs,type)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
rou=m/n;
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
        
        %%% if debug
        if xs~=-1
            FF=~AAz;
            Bu=Au(FF);
            LLB=Bu'*Bu;
            
            %Au=L(xk,zk)-L(xk+1,zk)
            Auu=Au;
            %Au-r;
            Aur=Auu-r0;
            %||Au||^2=||L(xk,zk)||^2-||L(xk+1,zk+1)||^2
            % Auun=Au(AA)
            Auun=Aur(AAz)'*Aur(AAz);
            Auunf=Aur(AAz);

            Auunff=Auunf'*Auunf;
            AuunfFF=Aur(FF);
            %AuunfFF(AuunfFF<0)=0;
            % AuunfFF=Au(FF)
            AuunfFFff=AuunfFF'*AuunfFF;
            %LMv=Ax_{k+1}-b-zk=(Au-r0+) 
            %LMv1=LMv(AA)=(Ax_{k+1}-b)(AA)=(Ax_{k+1}-b)
            LMv1=LMv(rkp0>ee);
            LMv1f=LMv1'*LMv1;
            %LMv2=LMv(FF)=(Ax_{k+1}-A_k)(FF)
            LMv2=LMv(rkp0<ee);
            LMv2f=LMv2'*LMv2;
            % if exposed point LMv1=L21 L22=Lmv2
            L21=rkp(rkp0>ee);
            L21f=L21'*L21;
            L22=rkp(rkp0<ee);
            L22(L22<0)=0;
            l22f=L22'*L22;
            LLMv1=LMv1f-L21f;
            LLMv1f=LLMv1'*LLMv1;
            LLMv2=LMv2-L22;
            LLMv2f=LLMv2'*LLMv2;
            
            sump=sum(rkp0>0);
            rks=(b-A*xs);
            sumpx=sum(rks>0);
      
            
            
            
            AAA=(rks>-1e-10);
            AII=A(AAA,:);
            ss=AII'*AII*xs-AII'*b(AAA);% valid xs true solution
            lll=min(eig(AII'*AII));% is positive then unqiue
            [xkkry,~]=krylov(A,b,xk,rkp);
            
            if  abs(LLA-LL)<diff
                fprintf('xs:%d,cs:%d,check ok\n',sumpx,sump);
            %else
            %    fprintf('xs:%d,cs:%d,check no\n',sumpx,sump);
            end
%             if sumpx==sump
            % 1：输入用于判断的参数 2：实际上收缩率 3：在充分大以后收缩率 4.最优点积极集 5 当前积极集 6，充分大后收缩量公式
            %7 解是否唯一 8输入的解是否为真 9-10 子空间方法是否缩短了距离
            fprintf('gradient|LL/LLA:%g,L11:%g,L12:%g,L21:%g,L22:%g,L211:%g,L222:%g,cs:%d,ac:%d,unqiue:%g,xs err:%g,dd1:%g,dd2:%g\n',...
                abs(LL-LLA),LLA,LLB,Auun-Auunff,LLB-AuunfFFff,LLMv1f,LLMv2f,sump,sumpx,lll,norm(ss),norm(xkkry-xs),norm(xk-xs));
 %            end
        end
        %%%
        
        %ssign=getBnS(eIter,Qn,fm,I);
        % if all great zeros mean same sign
        % if LLA>rou*LL || LLA<elta
        if abs(LLA-LL)<diff
            %newtonalgorithm
            countNW=countNW+1;
            % record begin newtron type iteratror
            if countNW ==1
                beginNW=countFM;
            end
            if type ==1
                [xk,~]=krylov(A,b,xk,rkp);
            else
                I=find(rkp>=ee);
                % 提取子矩阵判断是否正定
                AI=A(I,:);
                hk=AI\rkp(I);
                aa=piecewise(A,b,hk,xk);
                xk=xk+aa*hk;
            end
            rkp=b-A*xk;
            rk=rkp;
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
    % test
    if maxIter < countFM
        break;
    end
end
tf=etime(clock,t);
vk=sum(sign(rk));
%disp(['QQ:',num2str(var)]);
%disp(['hybrid6 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);