function [A,b,x0] = readData(type,m1,m2,n)
switch type
    case 1
        A1=abs(sprand(m1,n,0.1,1/100));
        A2=abs(sprand(m2,n,0.1,1/100));
        %A1=rand(m1,n) + 1 ;
        %A2=rand(m2,n) + 1;
        b1=rand(m1,1);
        b2=rand(m2,1);
        A=[A1;-A2];
        b=[b1;-b2];
        x0=ones(n,1);
   case 2
        A1=abs(rand(m1,n))/100;
        A2=abs(rand(m2,n))/100;
        %A1=rand(m1,n) + 1 ;
        %A2=rand(m2,n) + 1;
        b1=rand(m1,1);
        b2=rand(m2,1);
        A=[A1;-A2];
        b=[b1;-b2];
        x0=ones(n,1);

    case 3
        A = [1 3;2 4;-5 -6]; b = [5;6;-3]; x0 = [1;1];m1=2;m2=1;n=2;
    case 4
         load('testFbug')
%          m1 = 5;m2=5;n=7;
    case 5 
        load('testface')
        m1 = 300;m2=300;n=200;x0=ones(200,1);
    case 6 
        load('testalpha1')
        m1 = 1000;m2=1000;n=300;x0=ones(300,1);
    % pcg 
    case 7
%         A1=abs(sprand(m1,n,0.1,1/100));
%         A2=abs(sprand(m2,n,0.1,1/100));
        A1=rand(m1,n) + 1 ;
        A2=rand(m2,n) + 1;
        b1=rand(m1,1);
        b2=rand(m2,1);
        A=[A1;-A2];
        b=[b1;-b2];
%         [u,s,v] = svd(A);
%         s(1:n,1:n) = diag(rand(n,1));
%         A = u*s*v';
%         
%         L = ichol(A'*A);
%         A = L*L';
         x0=ones(n,1);
    otherwise
         load('test15')
         m1 = 500;m2=500;n=300;x0=ones(300,1);
end
%  hold on 
% clc
% clear
% load('testAbx0toldettoumaxitforconvergenceipg.mat');
