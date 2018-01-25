H=[4 5 -5;5 9 -5;-5 -5 7];
c=[2;1;-3];
lx=[-0.5;0;0];
n=3;
x0=-1*H\c;
z0=[0,0,0]';
x1=-1*H\c;
z1=[0,0,0]';
% while norm(x1-x0)+norm(z1-z0)>1e-5
% for i=0:5
% x0=x1;
% z0=z1;
% format rat
% x0=[0;0;0];
% z0=[-2;-1;3];
% [x1,z1,h]=semismooth(H,c,x0,z0,n)
% t=[x0;z0]./h
% t(t>-0.001)=-inf;
% t1=-1*max(t)
% [x0;z0]+t1*h

% HH=[
%        4              5             -5              1              0              0       
%        5              9             -5              0              1              0       
%       -5             -5              7              0              0              1       
%        0              0              0              1              0              0       
%        0              0              0              0              2/5              0       
%        0              0             -1              0              0              0  ];
%   [x0;z0]+ HH\[0,0,0,2,1,0]'
   
%disp(['x0:',num2str(x0'),' z0:',num2str(z0'),' x1:',num2str(x1'),' z1:',num2str(z1')]);
% end
% N=3;
% X = diag(rand(N,1));
% U = orth(rand(N,N));
% H = U' * X * U;
% x=rand(N,1)-0.5;
% c = -1*H*x;
format rat
% lx=quadprog(H,c,[],[],[],[],[],zeros(N,1))
lx(lx>-1e-8)=0
newtest(H,c,lx);