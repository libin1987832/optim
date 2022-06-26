
A=rand(5, 5);
charA = '11111';
dec2bin(5,5)==charA
iter = 1;
A=[5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
charA = '1111';
% B=[2 6 -1 -2;5 -1 2 3;-3 -4 1 10;5 -2 -3 8];
%[V,D]=eig(A,B)
for iter = 1 : 2^4-1
    F = dec2bin(iter, 5) == charA;
    AF = A(F, F);
    eig(A)
end
    