%%  BAS本程序目的：求出当A为对称（正定）矩阵，B为单位矩阵的特征值互补问题.
function [i] = spBas(A, B, x0, beta, etaMin, eps)
e=ones(n,1);      %生成n维的1向量
%Beta=10^(-5);
%varepsilon=10^(-5);%精度
%%  开始循环
k=0;
%Etamin=unifrnd (0,1);
etaMax=1/etaMin;

typeB = 1;
if size(B,2) == 1
    typeB = 2;
end
x = x0;
xp = x0;
dfp = x0;

while 1
    
    if k == 0 
        eta = rand(0,1) * etaMax;
    else
        s = x - xp;
        w = df - dfp;
        sw = s' * w;
        if sw <= 0
            eta = etaMax;
        else
            maxetas =max( etaMin, s' * s / sw );
            eta = min(etaMax, maxetas);
        end
    end
    Ax = A * x;
    Bx = computBx(B , x, typeB);
    xAx = x' * Ax;
    xBx = x' * Bx;
    invxBx = 1 / xBx;
    df = 2 * invxBx * ( xAx * invxBx * Bx - Ax);
    
    %L = x <= beta * df;
    d = -x;
    %etabeta = max(eta, beta);
    Fg = x >= eta * df && x > beta * df;
    d( Fg ) = -eta * df( Fg );
    
    if norm( d ) <= varepsilon   %判断精度终止
        break
    end
    %% 计算Deltak
    Ad = A * d;
    dAd = d' * Ad;
    xAd = x' * Ad;
    Bd = computBx(B , d, typeB);
    dBd = d' * Bd;
    xBd = x' * Bd;
    a0= xAd * xBx - xBd * xAx;
    a1= dAd * xBx - dBd * xAx;
    a2= dAd * xBd - dBd * xAd;
    b = sqrt( a1^2 -4 * a2 * a0 );
    root1 = 0.5 * ( -a1 - b ) / a2;
    root2 = 0.5 * ( -a1 + b ) / a2; %%二次函数的两个根
    deta1 = min (root1, root2);
    deta2 = max (root1, root2);
    detaf = @( deta ) (xAx + deta^2 * dAd +2 * deta * xAd) / (xBx + deta^2 * dBd + 2 * deta * xBd);
    %%Deltak赋值
    F0 = detaf(1) ;
    F1 = detaf(deta1);
    F2 = detaf(deta2);
    if Delta2>0&&Delta2<=1
        if Delta1>0&&Delta1<=1
            if F2<F0&&F2<F1
                Deltak=Delta2;
            elseif F1<F0&&F1<F2
                Deltak=Delta1;
            else
                Deltak=1;
            end
        elseif F2<F0
            Deltak=Delta2;
        else
            Deltak=1;
        end
    elseif Delta1>0&&Delta1<=1
        if F1<F0
            Deltak=Delta1;
        else
            Deltak=1;
        end
    else
        Deltak=1;
    end
    x0=x1;
    x3=x1+Deltak*dk;
    Muk=1/(e'*x3);
    x1=x3.* Muk;
    k=k+1;
    %% 计算Etak
    sk=x1-x0;
    f1x0=f1x1;
    L1=x1'*x1;
    L3=x1'*A*x1;
    f1x1=-(2/L1).*((L3/L1).*x1-A*x1); %f1表示函数F的导数
    wk=f1x1-f1x0;
    if sk'*wk<=0
        Etak= Etamax;
    else
        Etak= min (Etamax ,max(Etamin,(sk'*sk)/(sk'*wk)));
    end
end
%% 输出
toc;                       %计算所用时间
                        
s=toc;
lamda=(x1'*A*x1)/(x1'*x1) %特征值
b=A*x1-lamda*x1;
h=x1'*b;
i=k           %迭代次数

end

function Bx = computBx(B , x, type)
        if type == 1
            Bx = B * x;
        elseif type == 2
            Bx = B .* x;
        else
            Bx = x;
        end
end