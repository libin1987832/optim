% %°¸Àý1
% %
addpath('./IPG')
addpath('./exact')
addpath('./hybrid')
addpath('./subproblem')
clear 
clc
[A,b,x0] = readData(type);

options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
options.Display = 'off';
options.StepTolerance = 1e-13;
% [xk1, resvec, arvec,face1v,face2v, tf1]=fixedMatrix(A,b,x0,maxIterA*5,1e-15,options);
% [rpk1, normr1, xmin1, Ar, normKKT1 , face11, face21] = kktResidual(A, b, xk1 , [], 1);
% h=semilogy(beginp:maxIterA,arvec(beginp:maxIterA),'b+');
% h.LineStyle = '--';
% hold on k1
[rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, [0;57/61] , [], 1); 
 fprintf('& %s & %g & %g & %g & %g &%g /n','kk',normr0,xmin0,normKKT0,min(Ar0));

% maxIterA = 300;
tic
[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m1+m2)],b,[],[],[],[],zeros(n+m1+m2,1),Inf*ones(n+m1+m2,1),x0,options);
 tf0=toc;
 xk0 = x1(1:n);
 [rpk0, normr0, xmin0, Ar0, normKKT0 , faceX0, faceA0] = kktResidual(A, b, xk0 , [], 1); 
 fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ls',normr0,xmin0,normKKT0,min(Ar0),tf0);
return 
 % [xk1, resvec, arvec, face1vec, face2vec, tf1]=fixedMatrix(A,b,x0,maxIterA,1e-15,options);
% [rpk1, normr1, xmin1, Ar1, normKKT1 , faceX1, faceA1] = kktResidual(A, b, xk1 , [], 1); 
% %(A, b, x0, tol, nf, maxit, options, type)

% maxIterA = 200;
% [xk2,resvec2,arvec2,face1h,face2h,tf2] = hybridnnls(A,b,x0,1e-5,3, maxIterA, options, 'IPG');
% [rpk2, normr2, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, xk2 , [], 1); 
% %arvec2(arvec2 ==0 ) = mean(arvec2);
% h=semilogy(1:2:2*maxIterA,arvec2(1:2:2*maxIterA),'r+');
% hold on
% semilogy(2:2:2*maxIterA,arvec2(2:2:2*maxIterA),'bo')
%  fprintf('& %s & %g & %g & %g & %g &%g\n','Hybrid',normr2,xmin2,normKKT2,min(Ar2),tf2);%
 
 AI = A([1,5,7,10],[1,3,5,6,7]);
 bI = b([1,5,7,10]);
%  xs = AT\bI;
 xs = lsqminnorm(AI, bI);
 xss = zeros(7,1);
 xss([1,3,5,6,7]) = xs;
[rpks, normrs, xmins, Ars, normKKTs , faceXs, faceAs] = kktResidual(A, b, xss , [], 1); 
 fprintf('& %s & %g & %g & %g & %g \n','xs',normrs,xmins,normKKTs,min(Ars));%
 AI*x1([1,3,5,6,7])-bI
%  hold on 
% clc
% clear
% load('testAbx0toldettoumaxitforconvergenceipg.mat');
maxIterA = 400;
% [xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-13, 1e-8, 1-1e-10, 2*maxIterA,'IPG');
% [rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
%  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','IPG',normr3,xmin3,normKKT3,min(Ar3),tf3);
% for i = 20:30
%     tou = i/31;
% [xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-8, 1e-5, 20/31, maxIterA,'IPG');
% [rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
%  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','IPG',normr3,xmin3,full(normKKT3),min(Ar3),tf3);
% semilogy(1:maxIterA,arvec3(1:maxIterA),'b+');
 % end
% [xk3, resvec3, arvec3, faceXvec3, tf3]  = IPG(A, b, x0, 1e-13, 1e-8, 1-1e-10, 2*maxIterA,'ST');
% [rpk3, normr3, xmin3, Ar3, normKKT3 , faceX3, faceA3] = kktResidual(A, b, xk3, [], 1);
%  fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','ST',normr3,xmin3,normKKT3,min(Ar3),tf3);

%  h=semilogy(1:maxIterA,resvec3(1:maxIterA),'b+');

%[tf2 normKKT2 normr2
 %tf3 normKKT3 normr3]
% A'*A<0
%
% fprintf('& %s & %g & %g & %g & %g &%g \\\\\n','FM',normr1,xmin1,normKKT1,min(Ar1),tf1);
% fprintf('& %s & %g & %g & %g & %g &%g
% \\\\\n','Hybrid',normr2,xmin2,normKKT2,min(Ar2),tf2);%
      
%[norm(xk1-xk2)]
%[normKKT1 normKKT2;normr1 normr2]
%tf1
% face1 = face1(face1>0);
% face2 = face2(face2>0);
% face1f=fliplr(face1);
% face2f=fliplr(face2);
% face1f(1:10)
% face2f(1:10)

% A2(A2<0)=0.5;
% A1(A1<0)=1;
% A=[2,1;2,-1];
% b=[5;-3];
% [x,fk]=fsearchx(A,b,[1;1],0.000001,0.0001)
