load('lp_adlittle.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:100 150:50:350 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,300000,10^(-2),beta,lambda,'adlittle')


load('lp_agg.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:100 150:100:950 1000:1000:floor(m/1000)*1000 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,150000,10^(-2),beta,lambda,'agg')


load('lp_blend.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:100 150:50:350 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,3000000,10^(-3),beta,lambda,'blend')


load('lp_bandm.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:100 200:100:1500 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,3000000,10^(-2),beta,lambda,'bandm')



load('lp_brandy.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:100 150:50:300 400:100:1000 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,1,30000000,0.5*10^(-1),beta,lambda,'brandy')


load('lp_degen2.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 50:50:250 300:100:2400 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,30000,10^(-2),beta,lambda,'degen2')


load('lp_finnis.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 50:50:250 300:100:1000 1500:500:3000 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,30000000,5*10^(-2),beta,lambda,'finnis')


load('lp_recipe.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:200 250:50:550 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,300000,0.2*10^(-2),beta,lambda,'recipe')


load('lp_scorpion.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:40 50:50:450 500:100:1600 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,50000,0.5*10^(-2),beta,lambda,'scorpion')


load('lp_stocfor1.mat'); 
n = size(C,2);
m = size(C,1);

beta=[1 10:10:40 50:50:550 m];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(C,d,10,200000,10^(-1),beta,lambda,'stocfor1')


