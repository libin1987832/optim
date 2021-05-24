addpath('../dataInequality/');
addpath('../FM/')
clc
clear
Arr=zeros(2*10,2+2*10);%m,n,t,nf,iffind);
index=0;
nfiff=zeros(1,2*10);
% m = 1000;
% n = 100;
for ma = 1000:1000:2000
    for n = 100:100:1000
    index=index+1;
    for t=1:10
        rangeMax = 2;
        rangeMin = -2;
        A = 2 * rand(ma , n)-1;
        b = 2 * rand(ma , 1)-1;
        x00 = zeros(n , 1);
        [nf,iffind]=generateNf(A,b,x00,200);
        nfiff(1,2*t-1:2*t)=[nf,iffind]; 
    end
    Arr(index,:)=[ma,double(n),nfiff];
    end
end
Arr