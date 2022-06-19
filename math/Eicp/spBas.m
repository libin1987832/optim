%%  BAS������Ŀ�ģ������AΪ�Գƣ�����������BΪ��λ���������ֵ��������.
function [i] = spBas(A, B, x0, beta, etaMin, eps)
e=ones(n,1);      %����nά��1����
%Beta=10^(-5);
%varepsilon=10^(-5);%����
%%  ��ʼѭ��
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
    
    if norm( d ) <= varepsilon   %�жϾ�����ֹ
        break
    end
    %% ����Deltak
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
    root2 = 0.5 * ( -a1 + b ) / a2; %%���κ�����������
    deta1 = min (root1, root2);
    deta2 = max (root1, root2);
    detaf = @( deta ) (xAx + deta^2 * dAd +2 * deta * xAd) / (xBx + deta^2 * dBd + 2 * deta * xBd);
    %%Deltak��ֵ
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
    %% ����Etak
    sk=x1-x0;
    f1x0=f1x1;
    L1=x1'*x1;
    L3=x1'*A*x1;
    f1x1=-(2/L1).*((L3/L1).*x1-A*x1); %f1��ʾ����F�ĵ���
    wk=f1x1-f1x0;
    if sk'*wk<=0
        Etak= Etamax;
    else
        Etak= min (Etamax ,max(Etamin,(sk'*sk)/(sk'*wk)));
    end
end
%% ���
toc;                       %��������ʱ��
                        
s=toc;
lamda=(x1'*A*x1)/(x1'*x1) %����ֵ
b=A*x1-lamda*x1;
h=x1'*b;
i=k           %��������

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