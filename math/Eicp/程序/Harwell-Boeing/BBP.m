
function [x, F, iter, ninf] = BBP(A, b, F, p, ninf1, epsbbp, strategy) 
iter = 0; 
[ m, n ] = size(A);
M = [ A -ones( m , 1 ); ones( 1 , m ) 0 ];
h = [ b ; -1 ];
N = 1 : n;

while 1
    iter = iter + 1;
    T = setdiff( N , F );
    Fn = union( F,  n + 1  );
    MFF = M( Fn , Fn );
    MTF = M( T , Fn );
    hF = h( Fn );
    hT = h( T );
    z = zeros( n + 1 , 1 );
    v = zeros( n + 1 , 1 );
    z( Fn ) = - MFF \ hF;
    v( T ) = hT + MTF *  z( Fn );
    x = z( 1 : n , 1 );
    ninf = sum( z(F) < -epsbbp ) + sum( v(T) < -epsbbp  );
   
    if ninf < ninf1
         break;
    end
    if strategy == 0 && iter == p
        break;
    end
    if ninf == 0
        break;
    elseif strategy == 1
         U = union( F( z(F) < -epsbbp ), T( v( T ) < -epsbbp ) );
        r = min( U );
        if ismember(r, F)
            F = setdiff( F, r );
        else
            F = union( F, r );
        end
    else
        F = union( setdiff(F, F( z(F) < -epsbbp )), T( v( T ) <= epsbbp) );
    end
end
end