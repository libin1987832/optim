addpath('../');
addpath('../../LCP')
addpath('../../LCP/MPRGP')
clear
clc
[A,b,x0] = readData(1,500,500,300);
[m,n] = size(A);
% AI = [A,-speye(m)];
% M = AI' * AI;
% q = AI' * b;
% [w,z,retcode] = LCPSolve(M,q,1e-10,20);%Lemka(M,q,1e-10,20);
% [s,i,c]=exponent(A,b,1e-10);
% for i = 1:c
% [rpk2, normr2, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, s(:,i) , [], 1);
% fprintf('%d, %g, %g \n',i,normr2,normKKT2);
% end
% x = z(1:n);
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter');
options.Display = 'off';
%options.StepTolerance = 1e-13;

tic;[x1,f1,residual,exitflag,output,ff] = lsqlin([A,-eye(m)],b,...
    [],[],[],[],zeros(n+m,1),Inf*ones(n+m,1),x0,options);tf0=toc;xk0 = x1(1:n);
bim = b-A * x1(1:n) > -1e-3;
bin = x1(1:n) > 1e-1;
x = solvels(A,b,bim,bin);
[rpk2, normr2, xmin2, Ar2, normKKT3 , faceX, faceA] = kktResidual(A, b, x , [], 1);
[rpk2, normr1, xmin2, Ar2, normKKT2 , faceX2, faceA2] = kktResidual(A, b, x1(1:n) , [], 1);
%   1.127491279294015   1.127492552105187   0.000000000000001
%   0.000033134791754 (4)1.154013900031632   1.154017560454029
%   0.000000000000000   0.000226358996217(5) 
%   0.541046088233500   0.541049496018813   0.000000000000001
%   0.000041149271351 (6) 0.908655184900057   0.908657901827569
%   0.000000000000001   0.000114152651913(7)

if normr2<normr1
  save('test16.mat','A','b','x0')
[normr2 normr1 normKKT3 normKKT2 faceX2 faceX faceA2  faceA]
sum(bim)
sum(bin)
end
