m = 2000;
n = 50;
A = 2 * rand(m , n)-1;
    b = 2 * rand(m , 1)-1;
    x0 = ones(n , 1);
    maxIter = 900;
    nf = 3;
    str = ['D','C','R','P'];
    [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
     rkh=b-A*xkh;
     rkh(rkh<0)=0;
    dh=norm(rkh);
     gh=norm(A'*rkh);
    fprintf('han:%g %4.4f\n',gh,tfh);
    I = eye(n);
    xxx = [];
for i = 1:1000
u = rand(n,1);
u = u.* xkh*5;
dht = b-A*u;
dh = dht;
dh2 = -dht;
dh2(dh2<0) = 0;
dh(dh<0) = 0;
dh3 = b - A * xkh;
h=A'*dh;
d1 = -(u-xkh)'*(A'*A)*h;
d2 = h'*h;
d3 = d1-d2;
d4 = A'*dh2+A'*(dh3);
xxx = [xxx [d1;d2;d3;-d4'*h;-h'*A'*dh2;-h'*A'*(dh3)]];
end
sum(xxx(1,:) < xxx(2,:))