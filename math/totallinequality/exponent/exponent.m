function [solution, index ,count] = exponent(A,b,tol)
[m,n] = size(A);
xz = zeros(n,1);
AT = A';
solution = zeros(n,n);
count = 0;
index = zeros(m+n,n);
for im = 1:(2^m - 1)
    bim = logical(de2bi(im,m));
    for in = 1:(2^n - 1)
        bin = logical(de2bi(in,n));
        xz(:) = 0;
        [Atm,Atn] = size(AT(bin , bim));
        if Atn < Atm
            continue;
        end
        xz(bin) = (AT(bin , bim) * A(bim,bin))\(AT(bin , bim) * b(bim));
         r = b - A * xz;
         rb = r( ~bim );
         r(r < 0) = 0;
         Ar = AT * r;
%        if all(xz(in) > -tol) & all( AT * r < tol) & all(b( ~bim ) < A(~bim,:) * xz )
       if all(xz(bin) > -tol) &  all( rb < tol  ) & all(Ar < tol)
            count = count + 1;
            solution(:,count) = xz;
            index(:,count) = [bim,bin]'; 
        end
    end
end
