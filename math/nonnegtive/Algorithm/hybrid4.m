% Dax hybrid algorithm
function [xk,fk,xkArr,countFM,countNW,Q]=hybrid4(x0,A,b)
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

I=diag(ones(1,m));
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
while Ar>delt*rn && rn>delt
    r(r>0)=1;
    rk=diag(r);
    abs(norm(I-Q*Q') -1)
    abs(norm(I-Q*Q'*rk) -1)
    if 1
        countFM=countFM+1;
        %FM algorithm
        statFM=statFM+1;
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        xkArr=[xkArr;[xk',fk,0]];
        norm(r0)/norm(rk)
    else
        countNW=countNW+1;
        uIndex=0;
        %newtonalgorithm
        statSS=statSS+1;
        [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
        xkArr=[xkArr;[xk',fk,1]];
    end
    r=rk;
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
end 
disp(['AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM)]);