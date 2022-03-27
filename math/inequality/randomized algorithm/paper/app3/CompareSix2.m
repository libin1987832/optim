%%%%1
% x=zeros(1,640);
% x(1,60:170) = 0.5:0.5:6;
% x(1,230) = 9;
% x(1,280) = 7;
% x(1,380) = 2.1;
% x(1,450:550) = ones(1,11)*3.8;
% N=64000;
% x=spalloc(1,N,N);
% x(1,600:1700) = 5:(6-5)/1100:6;
% x(1,2300:3000) = 9;
% x(1,3800:4200) = 7;
% x(1,4300) = 2.1;
% x(1,4500:5500) = ones(1,1001)*3.8;
% N = length(x);%求取抽样点数
% t = (0:N-1);%显示实际时间
N=6000;
x=zeros(1,N);
x(1,60:170) = 5:(6-5)/110:6;
x(1,230:300) = 9;
x(1,380:420) = 7;
x(1,430) = 2.1;
x(1,450:550) = ones(1,101)*3.8;
N = length(x);%求取抽样点数
t = (0:N-1);%显示实际时间


% gausFilter = fspecial('gaussian',[1 3],3);
m=50;
 gausFilter = fspecial('gaussian',[1 m],3);
  L=zeros(N,N);
 for i = 1:m
   v1=  ones(1,N-i+1) * gausFilter(i);
   L1 = diag(v1,i-1);
   %L1=sparse(L1);
   L= L + L1;
% v1 = ones(1,N-1)*gausFilter(1);
% L1 = diag(v1,1);
% v2 = ones(1,N)*gausFilter(2);
% L2 = diag(v2);
% v3 = ones(1,N-1)*gausFilter(3);
% L3 = diag(v3,-1);
 end
% L = L1 + L2 + L3;
L=L+L';
 delt=0.15;
u=zeros(1,N);
uN=floor(N);
u(1:uN) = rand(1,uN)*0.9*delt;
x2 = L*x'+ u';
xu=12;
A=[-L;L;eye(N);-eye(N)];
b=[-x2-delt;x2-delt;zeros(N,1);-ones(N,1)*xu];
x0=zeros(N,1);

maxIter=200000;
iter = 30;
nf = 10;
tol = 1e-5;
xst=zeros(1,4);
t=clock;
[x1,arr]=project(L,x2,x0,iter,delt,0,xu,0);
xst(1)=etime(clock,t);
t=clock;
[x2,iter,error_k,iter_k,index_k] = GuassSeidelNE(A, b, x0,2 ,maxIter,tol,[],0);
xst(2)=etime(clock,t);
t=clock;
[x3,iter,error_k,iter_k,index_k] = simpleGuassSeidelNE(A, b, x0,2 ,maxIter,tol,[],0);
xst(3)=etime(clock,t);
t=clock;
[x4,iter,error_k,iter_k,index_k] = randGuassSeidelNE(A, b, x0,2 ,maxIter,tol,[],0);
xst(4)=etime(clock,t);

xs=[x1 x2 x3 x4];

SNRdB = @(s,n)( 10*log10(sum(s(:).^2)/sum((n(:)-s(:)).^2)) ); 
str=['P','C','U','R'];
for i = 1:4
r = b - A * xs(:,i);
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A'*r);
fprintf('& %sCD &%g &%g & %g & %g \\\\\n', str(i), SNRdB(x,xs(:,i)),r_GS, g_GS, xst(i));
end


% & PCD &17.9687 &0.182894 & 0.0986153 & 6.187 \\
% & CCD &17.0119 &0.524659 & 0.082917 & 5.673 \\
% & UCD &4.62 &1.60415 & 0.262694 & 10.364 \\
% & RCD &4.48563 &1.49949 & 0.231798 & 11.358 \\
% 
% & PCD &16.9358 &0.6605 & 0.47467 & 79.913 \\
% & CCD &16.6472 &0.675647 & 0.11328 & 57.519 \\
% & UCD &4.53101 &1.52782 & 0.287233 & 129.173 \\
% & RCD &4.30754 &1.96599 & 0.330704 & 142.153 \\

% 
% & PCD &16.907 &0.667337 & 0.513625 & 78.914 \\
% & CCD &15.5803 &1.28441 & 0.316829 & 29.142 \\
% & UCD &3.98662 &3.13809 & 0.826979 & 65.938 \\
% & RCD &3.96707 &3.78204 & 0.987099 & 72.053 \\