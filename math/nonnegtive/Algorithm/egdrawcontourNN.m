function egdrawcontourNN()
[X,Y] = meshgrid(-1:.05:0,0:.1:1);
Z=der(X,Y);
[C,h] = contour(X,Y,Z);
set(h,'ShowText','on','LevelStep',get(h,'LevelStep')*0.05)
colormap cool
hold on
% line([0,-1],[0,1])
line([-0.5,-0.5],[0,1])
line([-1,0],[0,2/3])
g=gradient(-0.5,1/3);
gn=newton(-0.5,1/3);
line([-0.5,g(1)-0.5],[1/3,g(2)+1/3]);
line([-0.5,gn(1)-0.5],[1/3,gn(2)+1/3]);
der(-0.5,1/3)
end
function Z=der(X,Y)
Z1=1-X-Y;
Z1(Z1<0)=0;
Z2=1+X+Y;
Z2(Z2<0)=0;
Z3=-0.5-X;
Z3(Z3<0)=0;
Z4=-2-2*X+3*Y;
Z4(Z4<0)=0;
Z=Z1.*Z1+Z2.*Z2+Z3.*Z3+Z4.*Z4;
Z=Z./2
end

function g=gradient(x,y)
A=[1,1;-1,-1;1,0;2,-3];
b=[1;1;-0.5;-2];
bA=b-A*[x;y];
bA(bA<0)=0;
g=A'*bA;
end

function gn=newton(x,y)
A=[1,1;-1,-1;1,0;2,-3];
b=[1;1;-0.5;-2];
bA=b-A*[x;y];
bA(bA<0)=0;
g=A'*bA;
gn=inv(A'*A)*g
end