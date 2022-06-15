function [d, lambda] = searchdir(M, Ax, xAx, xB, xBx, tao, D, x, method)
    if method == 1
        %tao��ֵ
      %%  ��dk,ֱ�ӵ��ú��� quadprog��(A�Գ�����ʱȡ��M=A.A�Գ�ʱȡM=A+ R*I)
       [d,fval,exitflag,output,lambda]=quadprog(M, Ax, [], [], xB, tao, -x, []);
       Md = M * d
       xMd = x' * Md;
       dMd = d * Md;
       dAx = d' * Ax;
       dBx = d' * xB';
       lambda = ( xMd + xAx + dMd + dAx )/( xBx + dBx);
    elseif method == 2
         Bx = xB';
         Dx = D .* x;
         maxBx = max(Bx(Bx>0));
         lBx = Bx/maxBx;
         lBxAx = lBx - Ax;
         grt0lBAx = lBxAx > Dx; 
         d = -x;
         d( grt0lBAx ) = lBxAx( grt0lBAx );
    end