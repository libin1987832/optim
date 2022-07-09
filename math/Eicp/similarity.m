function [flex, fed, fisd] = similarity(G,H)
 B = I;
M = I*3;
A = G;
for i = 1:2
    if i == 2
        A = H;
    end
epsx = 0;
epsxlambda = -1e-3;
x1 = ones(n,1)/n;
O=inv(B)*A;
R=max(abs(eig(O))); %AµÄÆ×°ë¾¶  
sigma0 = R + 0.01;
[xa, iter, error] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-8, 1000, 0);
nx = size(xa,2);
xd = []
for i = 1 : nx
    x = xa(:,i) ./ sum(xa(:,i));
lambda = (x' * A * x) / (x' * B * x)-20;
%xd = [xd [lambda; x]];
xd = [xd lambda];
%disp(['FQP_QUAD:lambda=' num2str(lambda) ', ninfx=' num2str(sum(x<epsx)) ',ninfy='  num2str(sum((A - lambda * B) * x < epsxlambda))])
end
 Glambda = sort(uniquetol(xd,1e-3),'descend');