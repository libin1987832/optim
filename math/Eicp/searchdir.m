function [d, lambda] = searchdir(M, Ax, xB, tao, D, x, method)
    if method == 1
        %tao��ֵ
      %%  ��dk,ֱ�ӵ��ú��� quadprog��(A�Գ�����ʱȡ��M=A.A�Գ�ʱȡM=A+ R*I)
       [d,fval,exitflag,output,lambda]=quadprog(M, Ax, [], [], xB, tao, -x, []);
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