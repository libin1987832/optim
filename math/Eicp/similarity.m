function [flex, fed, fisd] = similarity(G, H, GHtol, shiftLambda)
A = G;
maxL = 0;
for i = 1:2
    if i == 2
        A = H;
    end
    n = size(A,1);
    I = eye(n);
    B = I;
    M = I;
    x1 = ones(n,1)/n;
    R=max(abs(eig(A))); %AµÄÆ×°ë¾¶
    sigma0 = R + 0.01;
    [xa, ~, ~] = allsqp(A, B, M, x1, sigma0, 0.1, 1e-10, 1e-12, 1000, 0);
    nx = size(xa,2);
    lambda = zeros(nx,1);
    for j = 1 : nx
        x = xa(:,j) ./ sum(xa(:,j));
        lambda(j) = (x' * A * x) / (x' * B * x)-shiftLambda;
    end
     lambdaT = sort(uniquetol(lambda,0.01),'descend');
   %   sort(uniquetol(xd,0.01),'descend')'
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
    if abs(lambdaA(i,1)-lambdaA(i,2)) > GHtol
       if lex == 0
          lex = i;
       end
       simed = simed + 2^(-i);
       simisd = simisd + i^(-2);
    end
end
flex = 1 - 1 / lex;
fed = 1 - simed;
fisd = 1 - 6 * pi^(-2) * simisd;