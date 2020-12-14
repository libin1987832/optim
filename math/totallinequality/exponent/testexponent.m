addpath('../');
addpath('../../LCP')
addpath('../../LCP/MPRGP')
clear
clc
[A,b,x0] = readData(4,200,200,300);
[m,n] = size(A);
AI = [A,-speye(m)];
M = AI' * AI;
q = AI' * b;
[w,z,retcode,tf] = LCPSolve(M,q,1e-10,20);%Lemka(M,q,1e-10,20);

x = z(1:n);
[rpk1, normr1, xmin1, Ar1, normKKT1 , faceX, faceA] = kktResidual(A, b, x , [], 1);
fprintf('& %s & %g & %g & %g & %g & %g  \n','FM',normr1,xmin1,normKKT1,min(Ar1),tf); 
% [s,i,c] = exponent(A,b,1e-10);
% for i = 1:c
% [rpk2, normr2, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, s(:,i) , [], 1);
% fprintf('%d, %g, %g \n',i,normr2,normKKT2);
% end



% options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
% options.Display = 'off';
% %options.StepTolerance = 1e-13;
% 
% tic;[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
%     [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);tf0=toc;xk0 = x1(1:n);
% bim = b-A * x1(1:n) > -1e-3;
% bin = x1(1:n) > 1e-1;
% x = solvels(A,b,bim,bin);
% [rpk2, normr2, xmin2, Ar2, normKKT3 , faceX, faceA] = kktResidual(A, b, x , [], 1);
% [rpk2, normr1, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, x1(1:n) , [], 1);
% 
% 
% if normr2<normr1
%   save('test16.mat','A','b','x0')
% [normr2 normr1 normKKT3 normKKT2 faceX2 faceX faceA2  faceA]
% sum(bim)
% sum(bin)
% end
