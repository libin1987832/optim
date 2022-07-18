%用块主枢轴算法（BBP）求解凸二次规划问题
function [x, F, iter, ninf] = BBP2(Bx, A, n, a, xBx, maxItsub,  epsbbp) % epsilon2表示误差，e表示分量全为1的横向量，e0表示表示分量全为0的横向量，options为了不输出 quadprog中间的计算过程
iter = 0;  %迭代次数
M = [A -Bx -ones( n , 1 ); Bx' 0 0 ; ones( 1 , n ) 0 0]; % M矩阵 
h = [ - a * Bx;  -xBx; -1];
N = 1 : n+1;
F = 1 : n+1;
while 1
    iter = iter + 1;
    T = setdiff( N , F );
    Fn = union( F,  n + 2  );
    
    MFF = M( Fn , Fn );
    MTF = M( T , Fn );
    hF = h( Fn );
    hT = h( T );
    z = zeros( n + 2 , 1 );
    v = zeros( n + 2 , 1 );
    z( Fn ) = - MFF \ hF;
    v( T ) = hT + MTF *  z( Fn );
    x = z( 1 : n , 1 ); 
    ninf = sum( z(F) < -epsbbp ) + sum( v(T) < -epsbbp  );
    
    if iter == maxItsub
        break;
    end
    if ninf == 0
        break;
    else
        F = union( setdiff(F, F( z(F) < -epsbbp )), T( v( T ) <= epsbbp) );
    end
end
end




