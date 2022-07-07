%用块主枢轴算法（BBP）求解凸二次规划问题
%求解  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [x, F, iter, ninf, testwx] = BBP(A, b, F, maxIt, ninf1, eps, strategy, debug) % epsilon2表示误差，e表示分量全为1的横向量，e0表示表示分量全为0的横向量，options为了不输出 quadprog中间的计算过程
iter = 0;  %迭代次数
[ m, n ] = size(A);
M = [ A -ones( m , 1 ); ones( 1 , m ) 0 ]; % M矩阵
h = [ b ; -1 ];
N = 1 : n;
testwx = 0;
if debug 
    testwx = zeros( 2 * n, 1 );
end
while 1
    iter = iter + 1;
    T = setdiff( N , F );
    Fn = union( F,  n + 1  );
    MFF = M( Fn , Fn );
    MTF = M( T , Fn );
    hF = h( Fn );
    hT = h( T );
  %  z = zeros( n + 1 , 1 );
  %  v = zeros( n + 1 , 1 );
  z=sparse( n + 1 , 1) ;
   v=sparse( n + 1 , 1) ;
    z( Fn ) = - MFF \ hF;
    v( T ) = hT + MTF *  z( Fn );
    x = z( 1 : n , 1 );
    ninf = sum( z(F) < -eps ) + sum( v(T) < -eps  );
    if debug 
        testwx( 1 : n ) = x;
        testwx( n+1 : 2 * n ) = v( 1 : n );
    end
    if ninf < ninf1
         break;
    end
    if iter == maxIt
        break;
    end
    if ninf == 0
        break;
    elseif strategy == 1
         U = union( F( z(F) < -eps ), T( v( T ) < -eps ) );
        r = min( U );
        if ismember(r, F)
            F = setdiff( F, r );
        else
            F = union( F, r );
        end
    else
        F = union( setdiff(F, F( z(F) < -eps )), T( v( T ) <= eps) );
    end
end
end