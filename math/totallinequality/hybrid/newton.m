function [xs,rpk] = newton(A,b,n,rpk,x0)
display = false;
tol = 1e-2;

AA = (rpk > -tol);
RR = (x0 > tol);
% x0(~RR) = 0;
% subspace
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
options.Display = 'off';
options.StepTolerance = 1e-13;

AI = A(AA, RR);
bI = rpk( AA );
AA = (rpk>-tol);
type = 'lsqlin';
switch(type)
    case 'lsqrmx'
        % % subspace
        lsqrTol= 1e-5;
        maxIter =10;
        [u,flag,relres,iter,resvec,lsvec,out] = lsqrmx(AI,bI,lsqrTol,maxIter,[],[],zeros(size(AI,2),1),A(:,RR),b,x0(RR),AA,3);
        % u = lsqminnorm(AI, bI);
    case 'lsqlin'
        ns=size(AI,2);
        [x,f1,residual,exitflag,output,ff] = lsqlin(AI,b(AA),...
            [],[],[],[],zeros(ns,1),Inf*ones(ns,1),x0(RR),options);
        u= x - x0(RR);
    case 'md'
        u = AI\bI;
end

p = zeros(n, 1);
p(RR) = u;
% debug for test.m nf = 2
[~, xs, aranges, retcode] = arraySpiece(A, b, x0, p);
[rpk, normr1, ~, Ar1, ~, face21, face22] = kktResidual(A, b, xs , [], 1);
if display
    [~, normr0, ~, g, normKKT, face11, face12] = kktResidual(A, b, x0 , [], 1);
%    if p'*g > 0
 %       assert(p'*g < 0)
 %   end
    fprintf('newton(lsqin): normBF(%g,%g),gBF(%g,%g),xa(%d,%d,%d,%d) u(%g)\n',...
        normr0,normr1,norm(g),norm(Ar1),face11, face21,face12, face22,norm(u));
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