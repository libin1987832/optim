function [d, lambda] = searchdir(M, Ax, xB, tao, D, x, method)
    if method == 1
        %tao的值
      %%  求dk,直接调用函数 quadprog，(A对称正定时取，M=A.A对称时取M=A+ R*I)
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