% Dax hybrid algorithm and r+
function [xk,fk,xkArr,countFM,countNW,Q]=hybrid5(x0,A,b)
t=clock;
%compute hybrid uIter
[m,n]=size(A);
uIter=max(33,(m+n)/4);
%FM need a qr decompose
[Q,R]=qr(A);
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
statFM=0;
statSS=0;
xkArr=[];

countFM=0;
countNW=0;

I=diag(ones(1,m));
cont=1;
while cont==1
    if uIndex<uIter
        countFM=countFM+1;
        %FM algorithm
        statFM=statFM+1;
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
    else 
        
    countFM=countFM+1;
    %FM algorithm
    statFM=statFM+1;
    [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
    u=xk-x0;
    un=norm(u);
    cont=1;
    if norm(A*(u))<delt2*un || un<delt2
        rb=b-A*xk;
        rbI=rb;
        % Âú×ãÔ¼Êø
        rbI(rb<=delt3)=0;
        % ²»Âú×ãÔ¼Êø
        rbI(rb>delt3)=1;
        if sum(rbI)==0
            cont=0;
        else
            countNW=countNW+1;
            uIndex=0;
            %newtonalgorithm
            statSS=statSS+1;
            [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
            if lambe ==1
                rb2=b-A*(xk-x0);
                rbI2=rb2;
                % ²»Âú×ãÔ¼Êø
                rbI2(rb2<=delt3)=0;
                rbI2(rb2>delt3)=1;
                srba=max(rbI-rbI2);
                srbi=min(rbI-rbI2);
                if srba == srbi
                    cont=0;
                end
                xkArr=[xkArr;[xk',fk,1]];
            end
        end
        
    end
    r=rk;
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
end
tf=etime(clock,t);
disp(['hybrid5 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM),' cpu:',num2str(tf),' uIter:',num2str(uIter)]);