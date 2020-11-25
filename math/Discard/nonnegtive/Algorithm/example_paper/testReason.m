% 检测预测效果 返回1 积极集不固定 无法预测 返回2 估计值错误 real FM goundtrueth optimal estimate max
function [optimaln,real]=testReason(x00,A,b,n)
same=n;
[sm,sn]=size(A);
[Q,R]=qr(A);
Qn=Q(:,1:sn);
QQn=Qn*Qn';
rkn=b-A*x00;
Nk=rkn;
Nk(Nk>0)=1;
Nk(Nk<0)=0;
for i=1:n
    [xk,r0,rk,fk,fm,fr]=FM(x00,Q,R,A,b);
    x00=xk;
    rkn0=b-A*x00;
    Nk0=rkn0;
    Nk0(Nk0>0)=1;
    Nk0(Nk0<0)=0;
    % active set change break；
    if sum(xor(Nk0,Nk))~=0
        same=i-1;
        break;
    end
end
real=same;
if same > 0
    while(predict(QQn,A,b,x00,same,rkn,Nk)==-1&& same >0)
        same=same-1;
    end
    optimaln=same;
elseif same==0
    optimaln=0;
end

function r=predict(QQn,A,b,x0,n,rkn,Nk)
[sm,sn]=size(b);
I=diag(ones(sm,1));
Bn=I-QQn*diag(Nk);
rt=Bn^n*rkn;
rp=(I-n*QQn*diag(Nk))*rkn;
for i=1:sm
    if sign(rt(i))*sign(rp(i))<0
        r=-1;
        return
    end
end
r=1;
return