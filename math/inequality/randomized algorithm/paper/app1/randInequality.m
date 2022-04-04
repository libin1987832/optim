function [A,b,x0]=randInequality(m,n,rangeMax,rangeMin)
length = rangeMax - rangeMin;
A = length * rand(m , n) - rangeMin;
b = length * rand(m , 1) - rangeMin;
x0 = zeros(n,1);