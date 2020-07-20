%fh = @(x,y) [4,5]*[x;y]+0.5*[x,y]*[2 2;2 3]*[x;y];
x = linspace(-5,3);
m=size(x,2);
y = linspace(-4,1);
n=size(y,2);
[X,Y] = meshgrid(x,y);
z=zeros(m,n);
for i=1:m
    for j=1:n
    z(i,j) = [4,5]*[x(i,j);y(i,j)]+0.5*[x(i,j),y(i,j)]*[2 2;2 3]*[x(i,j);y(i,j)];
    end
end
contour(x,y,z)

%ezcontour(fh,[1,2,3]);