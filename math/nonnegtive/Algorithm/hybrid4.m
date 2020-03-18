% Dax hybrid algorithm and r+
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
am2=max(r);
ee=1e-15;% computer floating point arithmetic
delt=am*m*n*10*ee;
delt2=delt*100;
delt3=ee*max(b);
uIndex=0;
statFM=0;
statSS=0;
xkArr=[];

countFM=0;
countNW=0;

I=diag(ones(1,m));
%||A'(r)+||<=delt||(r)+|| ||(r)+||<=de
% while Ar>delt*rn && rn>delt
%     r(r>0)=1;
%     rk=diag(r);
%     abs(norm(I-Q*Q') -1)
%     abs(norm(I-Q*Q'*rk) -1)
while cont==1
        countFM=countFM+1;
        %FM algorithm
        statFM=statFM+1;
        [xk,r0,rk,fk,fm,fr]=FM(x0,Q,R,A,b);
        u=xk-x0;
        un=norm(u);
        cont=1;
    if norm(A(xk-x0))<delt2*un
        rb=b-xk;
        rbI=rb;
        % ������Լ��
        rbI(rb<=delt3)=0;
        rbI(rb>delt3)=1;
        countNW=countNW+1;
        uIndex=0;
        %newtonalgorithm
        statSS=statSS+1;
        [xk,rk,fk,f0,lambe]=ssqr(x0,A,b);
        rb2=b-xk;
        rbI2=rb2;
        % ������Լ��
        rbI2(rb2<=delt3)=0;
        rbI2(rb2>delt3)=1;
        srb=sum(rbI-rb2);
        if srb == 0
            cont=0;
        end
        xkArr=[xkArr;[xk',fk,1]];
    end
    r=rk;
    uIndex=uIndex+1;
    Ar=norm(A'*rk);
    rn=norm(rk);
    x0=xk;
end 
disp(['AT(b-A*x)+:',num2str(Ar),' fk:',num2str(fk),' ssqr:',num2str(statSS),' FM:',num2str(statFM)]);