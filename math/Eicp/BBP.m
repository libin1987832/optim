%�ÿ��������㷨��BBP�����͹���ι滮����
%���  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [x, F, iter]=BBP(A, b, F, maxIt, eps, strategy, debug) % epsilon2��ʾ��e��ʾ����ȫΪ1�ĺ�������e0��ʾ��ʾ����ȫΪ0�ĺ�������optionsΪ�˲���� quadprog�м�ļ������
iter = 0;  %��������
[ m, n ] = size(A);
M = [ A -ones( m , 1 ); ones( 1 , m ) 0 ]; % M����
h = [ b ; -1 ];
N = 1 : n;
% N=N';
F = N;

while 1
    T = setdiff( N , F );
    Fn = union( F,  n + 1  );
    MFF = M( Fn , Fn );
    MTF = M( T , Fn );
    hF = h( Fn );
    hT = h( T );
    z = zeros( n + 1 , 1 );
    v = zeros( n + 1 , 1 );
    z( Fn ) = - MFF \ hF;
    v( T ) = hT + MTF * zF;
    x = z( 1 : n , 1 );
    if iter > maxIt
        break;
    end
    if sum( z( F ) >= -eps ) == 0 && sum( v( T ) >= -eps) == 0 
        break
    elseif strategy == 1
        
    else
        F1 = setdiff( F , z < -eps );
        F2 = v( T ) <= eps;
        F = union(F1, F2);
    end
    iter = iter + 1;
end
end