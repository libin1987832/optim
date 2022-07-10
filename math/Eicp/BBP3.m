%用块主枢轴算法（BBP）求解凸二次规划问题
%求解  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [x, eta, F, iter, testwx] = BBP3(A, b, F, maxIt, eps, debug) % epsilon2表示误差，e表示分量全为1的横向量，e0表示表示分量全为0的横向量，options为了不输出 quadprog中间的计算过程
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
  %   z( Fn ) = - lsqminnorm(MFF,hF, 1e-30);
 %   z( Fn ) = -lsqr(MFF,hF,1e-3,100) ;
    %L = ichol(MFF,struct('michol','on'));
  %  L = ichol(MFF,struct('type','ict','droptol',10));
%    [x2,fl2,rr2,it2,rv2] = pcg(MFF,hF,1e-8,100,L,L');
%    z( Fn ) = -x2;
    v( T ) = hT + MTF *  z( Fn );
    x = z( 1 : n , 1 );
    eta = z(n+1);
    if debug 
        testwx( 1 : n ) = x;
        testwx( n+1 : 2 * n ) = v( 1 : n );
    end
    if size(F,2) > 0 && size(T,2) > 0 && min(z(F)) >= -eps && min(v(T)) >= -eps
        break;
    end
     if size(F,2)>0 && size(T,2) == 0 && min(z(F)) >= -eps 
        break;
     end
    if size(T,2)>0 && size(F,2) == 0 && min(v(T)) >= -eps
        break;
    end
    if iter == maxIt
        break;
    end
    F = union( setdiff(F, F( z(F) < -eps )), T( v( T ) <= eps) );
end
end