function [A,b,x0] = readData(type)
switch type
    case 1
        m1 = 1000;
        m2 = 1000; n = 200;
        A1=abs(sprand(m1,n,0.1,1/100))/100;
        A2=abs(sprand(m2,n,0.1,1/100))/100;
        %A1=rand(m1,n) + 1 ;
        %A2=rand(m2,n) + 1;
        b1=rand(m1,1);
        b2=rand(m2,1);
        A=[A1;-A2];
        b=[b1;-b2];
        x0=sparse(ones(n,1));
    case 2
         load('test')
         m1 = 5;m2=5;n=7;
    case 3
        A = [1 3;2 4;-5 -6]; b = [5;6;-3]; x0 = [1;1];m1=2;m2=1;n=2;
    otherwise
end
%  hold on 
% clc
% clear
% load('testAbx0toldettoumaxitforconvergenceipg.mat');
