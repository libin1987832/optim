x=-3:0.5:3;
y=-3:0.5:3;
[X,Y]=meshgrid(x,y);
z=0.5.*X.*X+X+2*Y;
% contour(X,Y,z)
%plot3(X,Y,z)
%surf(X,Y,z)
fsurf(@(x,y)05*x^2+x+2*y)