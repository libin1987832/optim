addpath('../IPG')
addpath('../exact')
addpath('../hybrid')
addpath('../subproblem')
addpath('../GNP')
addpath('../exponent')
addpath('../hybridfast')
addpath('../../inequality')
addpath('../../inequality/algorithmInequality')
addpath('../mprgp')
addpath('../class')
addpath('../')
clc
clear 
% 
[A,b,x0] = readData(1,200,200,200);

[m,n] =size(A);
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options = optimoptions('Algorithm','interior-point','TolX',1e-13)
options.Display = 'off';
% options.StepTolerance = 1e-13;
options.OptimalityTolerance = 1e-15;
options.ConstraintTolerance = 1e-20;
options.MaxIterations = 200;
% options.
param = parametern();
save('testMIL','A','b','x0','param','options');
% load('testMIL')
% options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
[m,n] =size(A);

example = 3;
if example == 3 || example > 10  
tic;[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
    [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);tf0=toc;xk0 = x1(1:n);
 [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, xk0 , [], 1); 
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ls',normr0,xmin0,normKKT0,min(Ar0),tf0);
 C=[A,-eye(m)]'*[A,-eye(m)];
%  0.5*x1'*C*x1-x1'*[A,-eye(m)]'*b
end

if example == 1 || example > 10
param.mprgp_a = 1/norm(A'*A, inf);
[xkn1, resvecn1, arvecn1,facen11v,facen12v, tfn1] = fixedMprgp(A,b,x0,param);
[rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xkn1 , [], 1);
fprintf('& %s & %g & %g & %g & %g & %g  \\\\\n','FMprgp',normr1,xmin1,normKKT1,min(Ar),tfn1);
end


if example == 2 || example > 10
maxIterA = 60000;
[xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-10, 1e-5, 0.8, maxIterA,'IPG');
[rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','IPG',normr3,full(xmin3),full(normKKT3),min(Ar3),tf3);
end 
if example == 3 || example > 10  
tic;
AE=[A,-eye(m)];
AEE = AE' * AE;
param.mprgp_a = 1/norm(AEE, inf)*2;
[x1,rpk] = MPRGPQ(AEE, AE'*b, [x0;ones(m,1)], param.mprgp_L, param.mprgp_a, param.mprgp_delta, param.mprgp_Ftol, 3000, b'*b);
tf0=toc;xk0 = x1(1:n);
 [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, xk0 , [], 1); 
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','mprgps',normr0,xmin0,normKKT0,min(Ar0),tf0);
 C=[A,-eye(m)]'*[A,-eye(m)];
%  0.5*x1'*C*x1-x1'*[A,-eye(m)]'*b
end

%  h=semilogy(1:maxIterA,resvec3(1:maxIterA),'b+');

%[tf2 normKKT2 normr2
 %tf3 normKKT3 normr3]
% A'*A<0
%

%% 
%  AI = A([1,5,7,10],[1,3,5,6,7]);
%  bI = b([1,5,7,10]);
% %  xs = AT\bI;
%  xs = lsqminnorm(AI, bI);
%  xss = zeros(7,1);
%  xss([1,3,5,6,7]) = xs;
% [rpks, normrs, xmins, Ars, normKKTs , faceXs, faceAs] = kktResidual(A, b, xss , [], 1); 
%  fprintf('& %s & %g & %g & %g & %g \n','xs',normrs,xmins,normKKTs,min(Ars));%
%  AI*x1([1,3,5,6,7])-bI
