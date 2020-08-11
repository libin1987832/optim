% Dax hybrid algorithm and r=B*r Nr==N nIter 预测的间隔 eIter 预测未来的多少步
function [xk,rk,countFM,countNW,beginNW,tf,vk]=contraction_i(x0,A,b,nIter,con,maxIter,xs,type)
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
        p1v=x0-xpn;
        p1=p1v'*p1v;
        p2v=xk-x0;
        p2=p2v'*p2v;
        roup=p2/p1;
        %%% if debug
        if xs~=-1
            AA=(rk>ee);
            IAA=eye(n)-pinv(A)*diag(AA)*A;
            lmax=max(eig(IAA));
            u1=p1v;
            u2=IAA*u1;
            rouI=(u2'*u2)/(u1'*u1);
            
           
            sump=sum(AA);
            rks=(b-A*xs);
            sumpx=sum(rks>0);
            AAA=(rks>-1e-10);
            FFF=(rks<1e-10);
            AII=A(AAA,:);
            NF=diag(FFF);
            ss=AII'*AII*xs-AII'*b(AAA);% valid xs true solution
            PAFA=pinv(A)*NF*A;
            p222=PAFA*p1v;
            rouI=max(abs(p222-p2));
            lll=max(eig(PAFA));% is positive then unqiue
            [xkkry,~]=krylov(A,b,xk,rkp);
            if roup<con
                fprintf('xs:%d,cs:%d,check ok\n',sumpx,sump);
            %else
            %    fprintf('xs:%d,cs:%d,check no\n',sumpx,sump);
            end
%             if sumpx==sump
            % 1：输入用于判断的参数 2：实际上收缩率 3：在充分大以后收缩率 4.最优点积极集 5 当前积极集 6，充分大后收缩量公式
            %7 解是否唯一 8输入的解是否为真 9-10 子空间方法是否缩短了距离
            fprintf('contraction|rou:%g,actural:%g,after expose:%g,xs:%d,cs:%d,ll:%g,unqiue:%g,xs err:%g,dd1:%g,  dd2:%g\n',...
                                    con,      roup,           rouI,sumpx,sump, lmax,       lll,norm(ss),norm(xkkry-xs),norm(xk-xs));
%             end
        end
        %%%
        
        %ssign=getBnS(eIter,Qn,fm,I);
        % if all great zeros mean same sign
        if roup<con
            %newtonalgorithm
            countNW=countNW+1;
            % record begin newtron type iteratror
            if countNW ==1
                beginNW=countFM;
            end
           % [xk,~]=krylov(A,b,xk,rkp);
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
%            if rou>0.6
%             nIter=(2^countNW)*nIter;
%             end
            %nIter=(2^countNW)*nIter;
            %             [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
            %             xkArr=[xkArr;[xk',fk,1]];
        end
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    xpn=x0;
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