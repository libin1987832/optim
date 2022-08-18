% load('creditdefault.mat');
% beta=[1 100:100:900 1000:1000:30000];
% lambda=[1 1.2 1.4 1.6 1.8];
% SKMtest(A,10,10000,10^(-2),beta,lambda,'Credit Card Default')

load('BreastCancerData.mat');
beta=[1 10:10:100 150:50:550 569];
lambda=[1 1.2 1.4 1.6 1.8];
SKMtest(A,10, 100000,0.0005,beta,lambda,'Wisconsin (Diagnostic) Breast Cancer')
