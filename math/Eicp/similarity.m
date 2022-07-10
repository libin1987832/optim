function [flex, fed, fisd] = similarity(G,H)
B = I;
M = I*3;
A = G;
maxL = 0;
for i = 1:2
    if i == 2
        A = H;
    end
    x1 = ones(n,1)/n;
    O=inv(B)*A;
    R=max(abs(eig(O))); %AµÄÆ×°ë¾¶
    sigma0 = R + 0.01;
    [xa, ~, ~] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-5, 1e-8, 1000, 0);
    nx = size(xa,2);
    lambda = zeros(nx,1);
    for j = 1 : nx
        x = xa(:,j) ./ sum(xa(:,j));
        lambda(j) = (x' * A * x) / (x' * B * x);
    end
     lambdaT = sort(uniquetol(lambda,1e-3),'descend');
     mT = size(lambdaT,1);
     if i == 1
        lambdaA = zeros(mT,2);
        lambdaA(:,1) = lambdaT;
        maxL = mT;
     elseif i == 2 && mT > maxL
         lambdaTT = lambdaA(:,1);
         lambdaA = zeros(mT,2);
         lambdaA(1:maxL,1) = lambdaTT;
         lambdaA(:,2) = lambdaT;
         maxL = mT;
     else
         lambdaA(1:mT,2) = lambdaT;
     end
end
simed = 0;
simisd = 0;
lex = 0;
for i = 1:maxL
    if abs(lambdaA(i,1)-lambdaA(i,2))>tol
       if lex == 0
          lex = i;
       end
       simed = simed + 2^(-i);
       simisd = simisd + i^(-2);
    end
end
flex = 1 - 1 / lex;
fed = 1 - simed;
fisd = 1-6 * pi * simisd;