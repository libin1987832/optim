% Dax hybrid algorithm and r=B*r Nr==N
function [xk,fk,xkArr,countFM,countNW,Q]=hybrid6(x0,A,b)
load('baseActive.mat');
t=clock;
%compute hybrid uIter
[m,n]=size(A);
uIter=max(33,(m+n)/4);
nIter=10;
%FM need a qr decompose
[Q,R]=qr(A);
Qn=Q(:,1:n);
QQ=Qn*Qn';
r=b-A*x0;
r(r<0)=0;
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
while Ar>delt*rn && rn>delt
    %     if uIndex<uIter
    if uIndex<nIter
        countFM=countFM+1;
        uIndex=uIndex+1;
        %FM algorithm
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
         % test
         if countFM<50
        rkt=b-A*xk;
        Nk2=rkt;
        Nk2(Nk2>0)=1;
        Nk2(Nk2<0)=0;
        reals=logical(Nk2'==active(countFM,:));
        realsum=sum(reals);
        if realsum~=m
            realsum
        end
         end
         % test 
    end
     % test 
    FAr=norm(A'*rk);
    Frn=norm(rk);
     % test 
    if mod(uIndex,10)==0 && uIndex >0
        % check if exposed face by r=B^n*r0
        uIndex=0;
        rkk=b-A*xk;
        Nk=rkk;
        Nk(Nk>0)=1;
        Nk(Nk<0)=0;
        Bn=I-nIter*QQ*diag(Nk);
        rkn=Bn*rk;
        % by corresponding compontent great zero
        signkn=sign(rkn.*rkk);
        % test 
        Nk2=rkk;
        Nk2(Nk2>0)=1;
        Nk2(Nk2<0)=0;
        reals=logical(Nk2'==active(base,:));
        realsum=sum(reals);
        % test
        ssign=sum(signkn);
        % if all great zeros mean same sign
        if ssign==m
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
    % test 
    if Ar>FAr
        [Ar,FAr]
    end
    % test 
    x0=xk;
end
tf=etime(clock,t);
disp(['hybrid6 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(countNW),' FM:',num2str(countFM),' cpu:',num2str(tf),' beginSS:',num2str(beginNW)]);