clear
clc
load('testpresc');
alpha1 = steplength * 1.02;%
maxa = max(alpha1,steplength);
mina = min(alpha1,steplength);
% knoty = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), [maxa,mina]);

xa = [ mina: (maxa-mina) / 100 : maxa];
ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
pxy={};
% pxy(1).X = [testalphax,testalphab];pxy(1).Y = knoty;
pxy(1).X = xa;pxy(1).Y = ya;
pxy(2).X = steplength;pxy(2).Y = funmin(A,b,x0,p,steplength);
pxy(3).X = alpha1;pxy(3).Y = funmin(A,b,x0,p,alpha1);
figure
hold on
p1 = arrayfun(@(a) plot(a.X,a.Y),pxy);
p1(1).Marker = 'o';p1(2).Marker = '.';p1(3).Marker = '*';%p1(4).Marker = 'x';
p1(1).LineStyle = 'none';
hold off


function f = funmin(A,b,x0,u,alpha)
xn = x0 + alpha*u;
xn( xn < 0) = 0;
f = b - A * xn;
f(f<0) = 0;
f = sqrt(f'*f);
end
