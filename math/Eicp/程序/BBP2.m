
function [x, F, iter, ninf] = BBP2(Bx, A, n, a, xBx, maxItsub,  epsbbp) 
iter = 0;  
M = [A -Bx -ones( n , 1 ); Bx' 0 0 ; ones( 1 , n ) 0 0]; 
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




