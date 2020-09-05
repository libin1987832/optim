% Dax hybrid algorithm and r=B*r Nr==N nIter 预测的间隔 eIter 预测未来的多少步
function [xk,rk,countFM,countNW,beginNW,tf,vk,rkArr]=scontraction(x0,A,b,nIter,con,maxIter)
t=clock;
rkArr=zeros(2*maxIter,1);
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
        p1v=x0-xpn;
        p1=p1v'*p1v;
        p2v=xk-x0;
        p2=p2v'*p2v;
        roup=p2/p1;
        
        % if all great zeros mean same sign
        if roup<con
            %newtonalgorithm
            countNW=countNW+1;
            % record begin newtron type iteratror
            if countNW ==1
                beginNW=countFM;
            end
           % [xk,~]=krylov(A,b,xk,rkp);
            [xk,~]=krylov(A,b,xk,rkp);
            rkp=b-A*xk;
            rk=rkp;
            rk(rk<0)=0;
        end
    end
    Ar=norm(A'*rk);
    rn=norm(rk);
    xpn=x0;
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