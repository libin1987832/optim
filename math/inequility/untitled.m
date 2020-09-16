%flops(0)           %start global flop count at 0
A=[1 2 3; 4 5 6];  
b=[7 8 9]';
x=A*b;             %do the operation
addflops(A*b) %do the counting
%flops              %display count so far