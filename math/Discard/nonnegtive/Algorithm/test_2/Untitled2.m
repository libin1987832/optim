% %fh = @(x,y) [4,5]*[x;y]+0.5*[x,y]*[2 2;2 3]*[x;y];
% x = linspace(-5,3);
% m=size(x,2);
% y = linspace(-4,1);
% n=size(y,2);
% [x,y] = meshgrid(x,y);
% z=zeros(m,n);
% b=[4;5];
% A=[2 2;2 3];
% for i=1:m
%     for j=1:n
%     z(i,j) = b'*[x(i,j);y(i,j)]+0.5*[x(i,j),y(i,j)]*A*[x(i,j);y(i,j)];
%     end
% end
% contour(x,y,z,[-4.3,-4,-4.4,-4.45,-4.49,-1,-4.1667],'ShowText','on')
% hold on 
% grid on
% x1=[-2;-1/3];
% A1=[2 0;0 15];
% b'*x1+0.5*x1'*A*x1
% b'*x1+0.5*x1'*A1*x1
%ezcontour(fh,[1,2,3]);
x1=Gauss(A,-b,[0;0],1e-5,1);
x2=Gauss(A,-b,[0;0],1e-5,2);
x3=Gauss(A,-b,[0;0],1e-5,3);
C=[x1 x2-x1 x3-x2]
D=A*[[0;0] x1 x2]+[b b b]
B=D(:,1:2)*inv(C(:,1:2));
B=[-5/3 -2;-2 -3];
B*C