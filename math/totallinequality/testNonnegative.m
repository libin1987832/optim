addpath('./IPG')
addpath('./exact')
addpath('./hybrid')
addpath('./subproblem')
addpath('./GNP')
addpath('./exponent')
addpath('./hybridfast')
addpath('../inequality')
addpath('../inequality/algorithmInequality')
clc
for i = 1:1
clear 

% 1:fixed 0:lsqin 2:hybrid lsqr 6:hybrd IPG 3: IPG 4: steep decent 5:Newton
% 7:GNP  8 hybridfast
example =1;
% 1 sparse matrix m1 m2 n 2 density matrix
[A,b,x0] = readData(2,1000,1000,500);
[m,n] =size(A);
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options = optimoptions('Algorithm','interior-point','TolX',1e-13)
options.Display = 'off';
% options.StepTolerance = 1e-13;
options.OptimalityTolerance = 1e-5;
% options.ConstraintTolerance = 1e-13;
options.MaxIterations = 600;
% options.
if example == 1 || example >100
maxIterA = 10;
[xk1, resvec, arvec,face1v,face2v, tf1] = fixedMatrix(A,b,x0,maxIterA,1e-15,options);
[rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xk1 , [], 1);
fprintf('& %s & %g & %g & %g & %g & %g  \n','FM1',normr1,xmin1,normKKT1,min(Ar),tf1); 
% [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, [0;57/61] , [], 1); 
%  fprintf('& %s & %g & %g & %g & %g &%g /n','FM',normr0,xmin0,normKKT0,min(Ar0));
%h=semilogy(beginp:maxIterA,arvec(beginp:maxIterA),'b+');
%h.LineStyle = '--';
 end
% 
if example == 1 || example >100
maxIterA = 10;
% options.OptimalityTolerance = 1e-10;
% options.Algorithm = 'active-set';
[xk1, resvec, arvec,face1v,face2v, tf1] = fixedMatrixlbfgs(A,b,x0,maxIterA,1e-15,options);
[rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xk1 , [], 1);
fprintf('& %s & %g & %g & %g & %g & %g  \n','FM2',normr1,xmin1,normKKT1,min(Ar),tf1); 
% [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, [0;57/61] , [], 1); 
%  fprintf('& %s & %g & %g & %g & %g &%g /n','FM',normr0,xmin0,normKKT0,min(Ar0));
%h=semilogy(beginp:maxIterA,arvec(beginp:maxIterA),'b+');
%h.LineStyle = '--';
end


%%
if example == 2 || example > 10  
tic;[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
    [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);tf0=toc;xk0 = x1(1:n);
 [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, xk0 , [], 1); 
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ls',normr0,xmin0,normKKT0,min(Ar0),tf0);
end
if example == 2 || example > 100  
 maxIterA = 100;
[xk2,resvec2,arvec2,face1h,face2h,tf2] = hybridnnls(A,b,x0,1e-5,3, maxIterA, options, 'ST');
[rpk2, normr2, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, xk2 , [], 1); 

% h=semilogy(1:2:2*maxIterA,resvec2(1:2:2*maxIterA),'r+');
%  hold on
%  semilogy(2:2:2*maxIterA,arvec2(2:2:2*maxIterA),'bo')
  fprintf('& %s & %g & %g & %g & %g &%g\n','HybridIsqr',normr2,xmin2,normKKT2,min(Ar2),tf2);%
end
if example == 6 || example > 100
    maxIterA = 10;
 [xk6,resvec6,arvec6,face6h,face6h,tf6] = hybridnnls(A,b,x0,1e-5,3, maxIterA, options, 'IPG');
[rpk6, normr6, xmin6, Ar6, normKKT6 , faceX6, faceA6] = kktResidual(A, b, xk6 , [], 1); 
% h=semilogy(1:2:2*maxIterA,resvec2(1:2:2*maxIterA),'r+');
%  hold on
%  semilogy(2:2:2*maxIterA,arvec2(2:2:2*maxIterA),'bo')
  fprintf('& %s & %g & %g & %g & %g &%g\n','HybridIPG',normr6,xmin6,normKKT6,min(Ar6),tf6);%
end



if example == 3 || example > 10
maxIterA = 60000;
[xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-10, 1e-2, 0.8, maxIterA,'IPG');
[rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','IPG',normr3,full(xmin3),full(normKKT3),min(Ar3),tf3);
end 
if example == 4 || example > 100
 [xk4, resvec4, arvec4, faceXvec4, tf4]  = IPG(A, b, x0, 1e-13, 1e-8, 1-1e-10, maxIterA,'ST');
 [rpk4, normr4, xmin4, Ar4, normKKT4 , faceX4, faceA4] = kktResidual(A, b, xk4, [], 1);
  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ST',normr4,full(xmin4),full(normKKT4),min(Ar4),tf4);
end
if example == 5 || example > 100
    maxIterA = 1;
 [xk5, resvec5, arvec5, faceXvec5, tf5]  = IPG(A, b, x0, 1e-13, 1e-8, 1-1e-10, maxIterA,'NT');
 [rpk5, normr5, xmin5, Ar5, normKKT5 , faceX5, faceA5] = kktResidual(A, b, xk5, [], 1);
  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','NT',normr5,full(xmin5),full(normKKT5),min(Ar5),tf5);
end
if example == 7 || example > 100
 M = 1e6;tol = 1e-4; delt = 1e-5;maxIterA = 1;
 [xkG, resvecG, arvecG, faceXvecG, tfG]  = GNP(A,b,x0,M,delt,tol,maxIterA);
 [rpkG, normrG, xminG, ArG, normKKTG , faceXG, faceAG] = kktResidual(A, b, xkG, [], 1);
  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','GNP',normrG,full(xminG),full(normKKTG),min(ArG),tfG);
end
if example == 8 || example > 10
    maxIterA =100;
 [xk8,resvec8,arvec8,face8h,face8h,tf8] = hybridfast(A,b,x0,1e-6,6, maxIterA);
[rpk8, normr8, xmin8, Ar8, normKKT8 , faceX8, faceA8] = kktResidual(A, b, xk8 , [], 1); 
% h=semilogy(1:2:2*maxIterA,resvec2(1:2:2*maxIterA),'r+');
%  hold on
%  semilogy(2:2:2*maxIterA,arvec2(2:2:2*maxIterA),'bo')
  fprintf('& %s & %g & %g & %g & %g &%g\n','Hybridfast',normr8,xmin8,normKKT8,min(Ar8),tf8);%
end
if example == 9 || example > 100
    maxIterA = 100;
 [xk9,resvec9,arvec9,face9h,face9h,tf9] = hybridprojnlss(A,b,x0,1e-6,6, maxIterA);
[rpk9, normr9, xmin9, Ar9, normKKT9 , faceX9, faceA9] = kktResidual(A, b, xk9 , [], 1); 
% h=semilogy(1:2:2*maxIterA,resvec2(1:2:2*maxIterA),'r+');
%  hold on
%  semilogy(2:2:2*maxIterA,arvec2(2:2:2*maxIterA),'bo')
  fprintf('& %s & %g & %g & %g & %g &%g\n','Hybridfastnnls',normr9,xmin9,normKKT9,min(Ar9),tf9);%
end
end
 % for i = 20:30
%     tou = i/31;
% [xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-8, 1e-5, 20/31, maxIterA,'IPG');
% [rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
%  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','IPG',normr3,xmin3,full(normKKT3),min(Ar3),tf3);
% semilogy(1:maxIterA,arvec3(1:maxIterA),'b+');
 % end

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
