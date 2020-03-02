% Dax hybrid algorithm
function [xk,fk,xkArr,countFM,countNW,Q]=hybrid1(x0,A,b)
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
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
uIndex=0;
statFM=0;
statSS=0;
xkArr=[];

countFM=0;
countNW=0;
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    if uIndex<uIter
        countFM=countFM+1;
        %FM algorithm
        statFM=statFM+1;
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
    else
        countNW=countNW+1;
        uIndex=0;
        %newtonalgorithm
        statSS=statSS+1;
        [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
        xkArr=[xkArr;[xk',fk,1]];
    end
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
    disp(['hybrid1 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM)]);
end 
tf=etime(clock,t);
disp(['hybrid1 m:',num2str(m),' n:',num2str(n),' AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM),' cpu:',num2str(tf)]);

    

    