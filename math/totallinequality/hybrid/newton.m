function [xs,rpk] = newton(A,b,n,rpk,x0)
display = true;
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
[rpk, normr, ~, g, normKKT, face1, face2] = kktResidual(A, b, x0 , rpk, 1);
if p'*g > 0
    assert(p'*g < 0)
end
% debug for test.m nf = 2
[alpha, aranges, retcode] = arraySpiece(A, b, x0, p);
xs = x0 + alpha * p;
xs(xs<0) = 0;
[rpk, normr, ~, g, normKKT, face1, face2] = kktResidual(A, b, x0 , rpk, 1);
if display
    fprintf('newton lsqin:p*g(%g),alpha(%g),minf(%g),KKT(%g)\n',p'*g,alpha,normr,normKKT);
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