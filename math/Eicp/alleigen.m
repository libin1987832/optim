
% A=rand(5, 5);
% charA = '11111';
% dec2bin(5,5)==charA
% iter = 1;
clc
clear
A=[5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
A = 0.5*(A'+A);
charA = '1111';
% B=[2 6 -1 -2;5 -1 2 3;-3 -4 1 10;5 -2 -3 8];
%[V,D]=eig(A,B)
alli = [];
allv = [];
alld = [];
epsx  = 0;
epsxlambda = -1e-10;
 for iter = 1 : 2^4-1
% for iter = 8:8
    F = dec2bin(iter, 4) == charA;
    AF = A(F, F);
    [v,d] = eig(AF);
    [m,n] = size(d);
    Fn = sum(F);
    for i = 1:n
        if sum(v(:,i) >= 0) == Fn
            v1 = zeros(4, 1);
            v1(F) = v(:,i) ./ sum(v(:,i));
            alli = [alli, iter];
            alld = [alld d(i,i)];
            allv = [allv v1];
            lambda = (v1' * A * v1) / (v1' * v1);
            disp(['lambda=' num2str(lambda) ', ninfx=' num2str(sum(v1<epsx)) ',ninfy='  num2str(sum(A *v1 - lambda * v1 < epsxlambda)) ',x(A-lb)x=' num2str(v1' * (A *v1 - lambda * v1))])
        end
    end
end
[alli;alld;allv]
min(abs(alld))
max(abs(alld))
    