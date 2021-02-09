addpath('./IPG');addpath('./exact');addpath('./hybrid');addpath('./subproblem');addpath('./GNP');addpath('./debug');
addpath('./exponent');addpath('./hybridfast');addpath('../inequality');addpath('../inequality/algorithmInequality');addpath('./mprgp');addpath('./class');
clc
clear 

for k = 1:1
% 1:fixed 0:lsqin 2: IPG 3:GNP 4:FMprgp
example = 11;
nn=[200,300,400,500];
% 1 sparse matrix m1 m2 n 2 density matrix
[A,b,x0] = readData(1,2000,2000,nn(k));[m,n] =size(A);



options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');options.Display = 'off';
options.OptimalityTolerance = 1e-15;options.ConstraintTolerance = 1e-10;options.MaxIterations = 70;
param = parametern();

 
if example == 0 || example > 10  
tic;[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
    [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);tf0=toc;xk0 = x1(1:n);
 [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, xk0 , [], 1); 
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ls',normr0,xmin0,normKKT0,min(Ar0),tf0);
end
if example == 1 || example >11
maxIterA = 200;
[xk1, resvec, arvec,face1v,face2v, tf1] = fixedMatrix(A,b,x0,maxIterA,1e-15,options);
[rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xk1 , [], 1);
fprintf('& %s & %g & %g & %g & %g & %g  \n','FM',normr1,xmin1,normKKT1,min(Ar),tf1); 
end
%%
%example = 3;
if example == 2 || example > 11
maxIterA = 5000;
[xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-6, 1e-5, 1, maxIterA,'IPG');
[rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','IPG',normr3,full(xmin3),full(normKKT3),min(Ar3),tf3);
end 
if example == 3 || example > 10
param.mprgp_a = 1/norm(A'*A, inf);
[xkn1, resvecn1, arvecn1,facen11v,facen12v, tfn1] = fixedMprgp(A,b,x0,param);
[rpk1, normr1, xmin1, Arn1, normKKT1 , face11, face21] = kktResidual(A, b, xkn1 , [], 1);
fprintf('& %s & %g & %g & %g & %g & %g  \\\\\n','FMprgp',normr1,xmin1,normKKT1,min(Arn1),tfn1);
end
end
