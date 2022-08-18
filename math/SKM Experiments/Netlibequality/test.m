
load('lp_agg.mat'); 
n = size(C,2);
m = size(C,1);

xinit = zeros(n,1);
skmeq(C,d,xinit,10,150000)



