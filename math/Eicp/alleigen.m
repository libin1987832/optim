
% A=rand(5, 5);
% charA = '11111';
% dec2bin(5,5)==charA
% iter = 1;
A=[5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
A = 0.5*(A'+A);
charA = '1111';
% B=[2 6 -1 -2;5 -1 2 3;-3 -4 1 10;5 -2 -3 8];
%[V,D]=eig(A,B)
alli = [];
allv = [];
alld = [];
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
            v1(F) = v(:,i);
            alli = [alli, iter];
            alld = [alld d(i,i)];
            allv = [allv v1];
            
        end
    end
end
[alli;alld;allv]
min(abs(alld))
max(abs(alld))
    