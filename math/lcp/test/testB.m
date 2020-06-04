%M=[2,1/2;1/2,1];
% M=[1/2,1;2,1/2];
% M=[1,1;2,1/2];
 M=[1/2,2;2,1/2];
q=[-1;-1];
x0=[0;2];
x1=[2/5;2/5];
x2=[2;0];
xr=rand(2,1)+1;
% x0=[1;0];
% x0=[1/3;2/3];
x0=[2;0];
[b,s]=checkB(M,q,x0)
[x,y] = meshgrid(-1:.05:2.2,-1:.05:2.2);
% x=2;
% y=0;
z = 0.25*x.^2+2.*x.*y+0.25*y.^2-x-y;
[C,h] = contour(x,y,z,50,'ShowText','on');
%  set(h,'ShowText','on','TextStep')
colormap cool
% 0.5*x1'*M*x1+q'*x1
% 0.5*x0'*M*x0+q'*x0
% 0.5*x2'*M*x2+q'*x2
% xr
% 0.5*xr'*M*xr+q'*xr