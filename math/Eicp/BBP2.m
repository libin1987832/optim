%�ÿ��������㷨��BBP�����͹���ι滮����
%���  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [x, F, iter, ninf, testwx] = BBP2(Ax, B, n, a, xAx, maxIt, ninf1, eps, strategy, debug) % epsilon2��ʾ��e��ʾ����ȫΪ1�ĺ�������e0��ʾ��ʾ����ȫΪ0�ĺ�������optionsΪ�˲���� quadprog�м�ļ������
iter = 0;  %��������
M = [ B -ones( n , 1 ) -Ax; ones( 1 , n ) 0 0; Ax' 0 0 ]; % M����
h = [ -Ax / a ; -1 ; -xAx];
N = 1 : n;
testwx = 0;
if debug 
    testwx = zeros( 2 * n + 1, 1 );
end
F = 1 : n;
while 1
    iter = iter + 1;
    T = setdiff( N , F );
    T = union(T, n+2);
    Fn = union( F,  n + 1  );
    
    MFF = M( Fn , Fn );
    MTF = M( T , Fn );
    hF = h( Fn );
    hT = h( T );
    z = zeros( n + 2 , 1 );
    v = zeros( n + 2 , 1 );
    z( Fn ) = - MFF \ hF;
    v( T ) = hT + MTF *  z( Fn );
    x = z( 1 : n , 1 ); 
    ninf = sum( z(F) < -eps ) + sum( v(T) < -eps  );
    if debug 
        testwx( 1 : n ) = x;
        testwx( n+1 : 2 * n ) = v( 1 : n );
        testwx( 2 * n + 1 : 2 * n + 1 ) = v( n + 2: n + 2 );
    end
    %if ninf < ninf1
    %     break;
   % end
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