H=[4 5 -5;5 9 -5;-5 -5 7];
c=[2;1;-3];
lx=[-0.5;0;0];
n=3;
x0=-1*H\c;
z0=[0,0,0]';
x1=-1*H\c;
z1=[0,0,0]';
% while norm(x1-x0)+norm(z1-z0)>1e-5
for i=0:5
x0=x1;
z0=z1;
[x1,z1,h]=semismooth(H,c,x0,z0,n);
disp(['x:',num2str(x1'),' z:',num2str(z1')]);
end
% N=3;
% X = diag(rand(N,1));
% U = orth(rand(N,N));
% H = U' * X * U;
% x=rand(N,1)-0.5;
% c = -1*H*x;
% lx=quadprog(H,c,[],[],[],[],[],zeros(N,1))
% lx(lx>-1e-8)=0
% newtest(H,c,lx);