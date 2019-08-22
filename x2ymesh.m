[X Y]=meshgrid(-5:0.1:5,-20:0.1:5)
Z=0.5*X.^2+X+Y;
[C,h]=contour(X,Y,Z);
colormap cool