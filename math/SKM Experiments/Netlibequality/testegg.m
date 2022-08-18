clc
clear
load('lp_agg.mat'); 
n = size(C,2);
m = size(C,1);
solution=SKM(A,b,xinit,error_bound,maxIter);