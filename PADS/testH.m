% H=[4 5 -5;5 9 -5;-5 -5 7];
% c=[2;1;-3];
% lx=[-0.5;0;0];
% n=3;
% N=3;
X = diag(rand(N,1));
U = orth(rand(N,N));
H = U' * X * U;
x=rand(N,1)-0.5;
c = -1*H*x;
lx=quadprog(H,c,[],[],[],[],[],zeros(N,1))
lx(lx>-1e-8)=0
newtest(H,c,lx);