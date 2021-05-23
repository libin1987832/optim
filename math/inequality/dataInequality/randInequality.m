function [A,b,x0]=randInequality(m,n,length,off,type)
if nargin==4 
    A = length * rand(m , n) - off;
    b = length * rand(m , 1) - off;
    x0 = ones(n,1);
else
    load(type);
end