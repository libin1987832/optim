function egdrawcontourNN()
[X,Y] = meshgrid(-1:.05:0,0:.1:1);
Z=der(X,Y);
% mesh(X,Y,Z)
[C,h] = contour(X,Y,Z);
set(h,'ShowText','on','LevelStep',get(h,'LevelStep')*0.05)
colormap cool
hold on
line([0,-1],[0,1])
line([-0.5,-0.5],[0,1])
% line([-1,0],[0,2/3])
line([-1/3,-5/6],[0,1])
g=gradient(-0.5,1/3);
gf=FM(-0.5,1/3);
gn=newton(-0.5,1/3);
H=line([-0.5,g(1)-0.5],[1/3,g(2)+1/3]);
 set(H,'Color','red')
H=line([-0.5,gf(1)-0.5],[1/3,gf(2)+1/3]);
 set(H,'Color','black')
line([-0.5,gn(1)-0.5],[1/3,gn(2)+1/3]);
der(-0.5,1/3)
end
function Z=der(X,Y)
Z1=1-X-Y;
Z1(Z1<0)=0;
Z2=1+X+Y;
Z2(Z2<0)=0;
Z3=0.5+X;
Z3(Z3<0)=0;
% Z4=-2-2*X+3*Y;
Z4=2+6*X+3*Y;
Z4(Z4<0)=0;
Z=Z1.*Z1+Z2.*Z2+Z3.*Z3+Z4.*Z4;
Z=Z./2
end

function g=gradient(x,y)
% A=[1,1;-1,-1;1,0;2,-3];
% b=[1;1;-0.5;-2];
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
bA=b-A*[x;y];
bA(bA<0)=0;
g=A'*bA;
end

function gf=FM(x,y)
% A=[1,1;-1,-1;1,0;2,-3];
% b=[1;1;-0.5;-2];
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
bA=b-A*[x;y];
bA(bA<0)=0;
g=A'*bA;
gf=inv(A'*A)*g
end

function gn=newton(x,y)
% A=[1,1;-1,-1;1,0;2,-3];
% b=[1;1;-0.5;-2];
A=[1,1;-1,-1;-1,0;-6,-3];
b=[1;1;0.5;2];
bA=b-A*[x;y];
NA=(bA>0);
AI=A(NA,:);
rI=b(NA)-AI*[x;y];
g=AI'*rI;
gn=(AI'*AI)\g
end