% m = 2000;
% n = 50;
% A = 2 * rand(m , n)-1;
%     b = 2 * rand(m , 1)-1;
A=[1 2;3 4];
b=[3;4];
m = 2;
n = 2;
x0 = ones(n , 1);
    maxIter = 900;
    nf = 3;
    str = ['D','C','R','P'];
   % [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
   xkh = [-2;2.5]; 
   rkh=b-A*xkh;
     rkh(rkh<0)=0;
    dhn=norm(rkh);
     gh=norm(A'*rkh);
     I = eye(n);
    xxx = [];
for i = 1:1000
u = 2*rand(n,1)-1;
y = 2*rand(m,1)-1;
y1 = y;
y2 = -y;
y1(y1<0) = 0;
y2(y2<0) = 0;
u = u.*xkh;
dht = b-A*u;
dh = dht;
dh2 = -dht;
dh2(dh2<0) = 0;
dh(dh<0) = 0;
dh3 = b - A * xkh;
dh5 = -dh3;
dh5(dh5<0) = 0;
h=A'*dh;
d1 = -(u-xkh)'*(A'*A)*h;
d2 = h'*h;
d3 = d1-d2;
d4 = A'*dh2+A'*(dh3);
%xxx = [xxx [d1;d2;d3;-d4'*h;-h'*A'*dh2;-h'*A'*(dh3);h'*A'*dh5]];
xxx = [xxx [d1;d2;d3;-d4'*h;h'*A'*dh5;y1'*A*A'*y2]];
end
xxx
sum(xxx(1,:) < xxx(2,:))
  fprintf('han:%g %4.4f %g\n',gh,tfh,dhn);
 
   
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