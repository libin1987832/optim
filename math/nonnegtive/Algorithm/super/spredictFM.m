% Dax hybrid algorithm and r=B*r Nr==N nIter 预测的间隔 eIter 预测未来的多少步
function [xk,rk,countFM,countNW,beginNW,tf,vk,rkArr]=spredictFM(x0,A,b,nIter,eIter,maxIter,xs,type)
% test baseActive.mat
% load('baseActive.mat');
t=clock;
rkArr=zeros(2*maxIter,1);
%compute hybrid uIter
[m,n]=size(A);
rou=m/n;
%FM need a qr decompose
[Q,R]=qr(A);
Qn=Q(:,1:n);
rkp=b-A*x0;
r=rkp;
r(r<0)=0;
fk=0.5*(r'*r);
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
while Ar>delt*rn && rn>delt
    %     if uIndex<uIter
    if uIndex<nIter
        countFM=countFM+1;
        uIndex=uIndex+1;
        %FM algorithm
        [xk,rkp]=FixedM(x0,Q,R,A,b,rkp);
        rk=rkp;
        rk(rk<0)=0;
    end
    
    if mod(uIndex,nIter)==0 && uIndex >0
        % check if exposed face by r=B^n*r0
        uIndex=0;
        %QQ=eIter*(Qn*Qn');
        AA=find(rk>ee);
        qrkn=Qn(AA,:)'*rk(AA);
        qrkn=Qn*qrkn;
        rkn=rkp-eIter*qrkn;
        ssign=sum(~xor(rk>ee,rkn>ee));
        if xs~=-1
        %%% if debug
        rkpN=rkp;
        xkN=xk;
        [xkN,rkpN]=FixedM(xkN,Q,R,A,b,rkpN);
        rknNN=rkp-qrkn;
        % predict formula is OK by one step for test Qn*Qn'
        nn=norm(rknNN-rkpN);
        rkpN=rkp;
        xkN=xk;
        % predict n step
        for i=1:eIter
            [xkN,rkpN]=FixedM(xkN,Q,R,A,b,rkpN);
        end
        sump=sum(rkn>0);
        sumpp=sum(rkpN>0);
        rks=b-A*xs;
        sumpx=sum(rks>0);
        AAA=(rks>-1e-10);
        AII=A(AAA,:);
        ss=AII'*AII*xs-AII'*b(AAA);
        lll=min(eig(AII'*AII));
        [xkkry,~]=krylov(A,b,xk,rkp);
         if  ssign==m ||ssign>m*0.99
                fprintf('xs:%d,cs:%d,check ok\n',sumpx,sump);
         %else
         %       fprintf('xs:%d,cs:%d,check no\n',sumpx,sump);
         end
%        if sumpx==sump
        % formula <1e-10 OK,predict active set ,FM true active set ,the
        % active set of solution ,active set of current point ,the
        % ununiqueness of solution,
        fprintf('predict|formula:%g,predict:%d,FM:%d,xs:%d,ss:%d,ll:%g,jj:%g,dd1:%g,dd2:%g\n', nn,sump,sumpp,sumpx,ssign,lll,norm(ss),norm(xkkry-xs),norm(xk-xs));
        %%%
 %       end
        end
        %ssign=getBnS(eIter,Qn,fm,I);
        % if all great zeros mean same sign
        if ssign==m ||ssign>m*0.99
            %newtonalgorithm
            countNW=countNW+1;
            % record begin newtron type iteratror
            if countNW ==1
                beginNW=countFM;
            end
            %[xk,~]=krylov(A,b,xk,rkp);
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
          %  nIter=countNW*nIter;
%           if rou>0.6
%              nIter=(2^countNW)*nIter;
%           end
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