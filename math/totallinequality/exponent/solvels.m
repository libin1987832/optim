function xz = solvels(A,b,bim,bin)
[ m , n ] = size(A);
xz = zeros(n,1);
AT = A';
xz(bin) = (AT(bin , bim) * A(bim,bin))\(AT(bin , bim) * b(bim));
% r = b - A * xz;
% rb = r( ~bim );
% r(r < 0) = 0;
% Ar = AT * r;
% %        if all(xz(in) > -tol) & all( AT * r < tol) & all(b( ~bim ) < A(~bim,:) * xz )
% if all(xz(bin) > -tol) &  all( rb < tol  ) & all(Ar < tol)
%     count = count + 1;
%     solution(:,count) = xz;
%     index(:,count) = [bim,bin]';
%     %       break;
% end
