function [xs,rpk] = newton(A,b,n,rpk,x0)
display = false;
tol = 1e-10;
 
AA = (rpk > -tol);
RR = (x0 > tol);
% x0(~RR) = 0;
% subspace
AI = A(AA, RR);
bI = rpk( AA );
 u = lsqminnorm(AI, bI);
% u = (AI+)\bI;
p = zeros(n, 1);
p(RR) = u;
% debug for test.m nf = 2
[~, xs, aranges, retcode] = arraySpiece(A, b, x0, p);
[rpk, normr1, ~, Ar1, normKKT, face1, face2] = kktResidual(A, b, xs , rpk, 1);
if display
    [rpk, normr, ~, g, normKKT, face1, face2] = kktResidual(A, b, x0 , rpk, 1);
    if p'*g > 0
        assert(p'*g < 0)
    end
    fprintf('newton(lsqin): normB(%g),normF(%g),gB(%g),gF(%g),xa(%d,%d),ba(%d,%d)',...
        normr0,normr1,Ar0,Ar1,face11, face21,face12, face22);
end
% rpk = b - A * xs;
% xa = [0:alpha/100:alpha*1.4];
% ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
%  pxy={};
% % pxy(1).X = knot;
% % pxy(1).Y = knoty;
% pxy(1).X = xa;
% pxy(1).Y = ya;
% figure
% p1 = arrayfun(@(a) plot(a.X,a.Y),pxy);
%hold on