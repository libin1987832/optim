% Dax hybrid algorithm and r=B*r Nr==N nIter 预测的间隔 eIter 预测未来的多少步 
function [xk,fk,xkArr,countFM,countNW,Q]=hybrid6(x0,A,b,nIter,eIter,maxIter,varargin)
% test baseActive.mat
% load('baseActive.mat');
var=size(varargin,2);
t=clock;
%compute hybrid uIter
[m,n]=size(A);
uIter=max(33,(m+n)/4);
% nIter=10;
%FM need a qr decompose
[Q,R]=qr(A);
Qn=Q(:,1:n);
if var==0 
 QQ=eIter*(Qn*Qn');
elseif var==2
    tmpq=3*ones(m,m);
end
r=b-A*x0;
r(r<0)=0;
fk=0.5*(r'*r);
%condition for terminate
Ar=norm(A'*r);
rn=norm(r);
am=max(max(A));
am2=max(r);
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
delt2=delt*100;
delt3=ee*max(b)*10^8;
uIndex=0;
xkArr=[];

countFM=0;
countNW=0;
beginNW=0;
I=diag(ones(1,m));
xk=x0;
while Ar>delt*rn && rn>delt
    %     if uIndex<uIter
    if uIndex<nIter
        countFM=countFM+1;
        uIndex=uIndex+1;
        %FM algorithm
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
    end

    if mod(uIndex,nIter)==0 && uIndex >0
        % check if exposed face by r=B^n*r0
        uIndex=0;

            if var==0
                ssign=getBn(QQ,fm,I);
            elseif var==1
                ssign=getBn2(eIter,Qn,fm,I);
            else
                [ssign,tmpq]=getBn2(eIter,Qn,fm,I,tmpq);    
            end
        % if all great zeros mean same sign
        if ssign==m ||ssign>m*0.99
            %newtonalgorithm
            countNW=countNW+1;
            % record begin newtron type iteratror
            if countNW ==1
                beginNW=countFM;
            end   
            [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
            xkArr=[xkArr;[xk',fk,1]];
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
disp(['hybrid6 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);