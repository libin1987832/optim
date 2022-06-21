%%  BAS������Ŀ�ģ������AΪ�Գƣ�����������BΪ��λ���������ֵ��������.
function [x, iter, error] = spBas(A, B, x0, beta, etaMin, eps, maxIt, debug)
iter = 0;
%Etamin=unifrnd (0,1);
etaMax=1/etaMin;
error = 0;
if debug
    error = zeros(1, maxIt + 1);
end
typeB = 1;
if size(B,2) == 1
    typeB = 2;
end
x = x0;
xp = x0;
dfp = x0;

while 1
    Ax = A * x;
    Bx = computBx(B , x, typeB);
    xAx = x' * Ax;
    xBx = x' * Bx;
    if debug
        error( iter + 1 ) = - xAx/xBx;
    end
    invxBx = 1 / xBx;
    df = - 2 * invxBx * ( xAx * invxBx * Bx - Ax);
    
     if iter == 0 
        eta = rand(1) * etaMax;
    else
        s = x - xp;
        w = df - dfp;
        sw = s' * w;
        if sw <= 0
            eta = etaMax;
        else
            maxetas = max( etaMin, s' * s / sw );
            eta = min(etaMax, maxetas);
        end
    end
    
    %L = x <= beta * df;
    d = -x;
    %etabeta = max(eta, beta);
    Fg = (x >= (eta * df)) & (x > (beta * df));
    d( Fg ) = -eta * df( Fg );
    
    if norm( d ) <= eps || iter > maxIt  %�жϾ�����ֹ
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
    detaf1 = detaf(1) ;
    detafdeta1 = detaf(deta1);
   % detafdeta2 = detaf(deta2);
    if deta2 < 0 || deta1 >1 
        deta = 1;
    elseif detaf1 < detafdeta1
        deta = 1;
    else
        deta = deta1;
    end
    xp = x;
    dfp = df;
    x = x + deta * d;
    Muk = 1 / sum(x);
    x = x .* Muk;
    iter = iter + 1;
end

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