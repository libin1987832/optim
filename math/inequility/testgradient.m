m = 50;
n = 10;
A = 2 * rand(m , n)-1;
% [q,r] =qr(A);
% D = diag(rand(n,1));
% A = D*q;
b = 2 * rand(m , 1)-1;
ATA=A'*A;
r=1/(max(eig(ATA))+0.0001);
min(eig(ATA))
x0 = ones(n , 1); 
[xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,100);
rkh = b- A*xkh;
rkh(rkh<0) = 0;
fxh = 0.5*rkh'*rkh;
arr = zeros(3,10);
% arr(1,:) = xkh*ones(1,10);
x = x0;
for i = 1 :30
bxt = b-A*x;
bx = bxt;
bx(bx<0) = 0;
fkhp = 0.5*bx'*bx;
d = A'*bx;
ddd = (r*A*A'-eye(m))*bx;
ffkhp = 0.5*ddd'*ddd;
 
xk = x + r*d;
rkhnt = b - A*xk;
rkhn = rkhnt;
rkkhntt = -rkhnt;
rkkhntt(rkkhntt<0) = 0;
bxn = -bxt;
bxn(bxn<0) =0;
rr = -rkhnt-bxn;
fff= 0.5*rr'*rr;
rkhn(rkhn<0) = 0;
fkhn = 0.5*rkhn'*rkhn;
x = xk;
arr(1,i)=fkhp;
arr(2,i)=fkhn;
arr(3,i)=(fkhn-fxh)/(fkhp-fxh);
end
arr

% x=x0;
% Ai=pinv(ATA);
% for i = 1 :10
% bx = b-A*x;
% bx(bx<0) = 0;
% d = A'*bx;
% xk = x + Ai*d;
% x = xk;
% arr(2,i)=norm(d);
% end
% arr


%-------------------
% A=[1 2;3 4];
% b=[3;4];
% m = 2;
% n = 2;
% x0 = ones(n , 1);
%     maxIter = 900;
%     nf = 3;
%     str = ['D','C','R','P'];
%    % [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
%    xkh = [-2;2.5]; 
%    rkh=b-A*xkh;
%      rkh(rkh<0)=0;
%     dhn=norm(rkh);
%      gh=norm(A'*rkh);
%      I = eye(n);
%     xxx = [];
% for i = 1:1000
% u = 2*rand(n,1)-1;
% y = 2*rand(m,1)-1;
% y1 = y;
% y2 = -y;
% y1(y1<0) = 0;
% y2(y2<0) = 0;
% u = u.*xkh;
% dht = b-A*u;
% dh = dht;
% dh2 = -dht;
% dh2(dh2<0) = 0;
% dh(dh<0) = 0;
% dh3 = b - A * xkh;
% dh5 = -dh3;
% dh5(dh5<0) = 0;
% h=A'*dh;
% d1 = -(u-xkh)'*(A'*A)*h;
% d2 = h'*h;
% d3 = d1-d2;
% d4 = A'*dh2+A'*(dh3);
% %xxx = [xxx [d1;d2;d3;-d4'*h;-h'*A'*dh2;-h'*A'*(dh3);h'*A'*dh5]];
% xxx = [xxx [d1;d2;d3;-d4'*h;h'*A'*dh5;y1'*A*A'*y2]];
% end
% xxx
% sum(xxx(1,:) < xxx(2,:))
%   fprintf('han:%g %4.4f %g\n',gh,tfh,dhn);
 
   
%  m = 10;
%  n = 2;
%  A = 2 * rand(m , n)-1;
%  AA = A*A';
%  x = (2*rand(m,1)-1)*100;
%  x1 = x;
%  x2 = -x;
%  y = (2*rand(m,1)-1);
%  y1 = y;
%  y2 = -y;
%  x1(x1<0)=0;
%  x2(x2<0)=0;
%  y1(y1<0)=0;
%  y2(y2<0)=0;
%  xy1 = [y1-x1]';
%  xy2 = [x2-y2];
% 
%  xy1*AA*xy2