clear
clc
debug = 1;
%% 产生问题矩阵


% 二维矩阵
%  A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
A = [2,3;3,-3];b=[2;3];x0=[1;1];
% 不一致情况下的正解
x_exact=[1;0];
% 一致情况下的正解
% x_exact=[0;0];

%% Kaczmarz
t=clock;
[x_Kac, iter_Kac, error_Kac, xA_Kac, index_Kac] = randomizedKaczmarzNE(A, b, x0, maxit_Rand,tol,x_exact,debug);
tf_Kac=etime(clock,t);
r = b - A * x_Kac;
r(r<0) = 0;
r_Kac = norm(r);
g_Kac = norm(A'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'Kaczmarz', r_Kac, g_Kac, iter_Kac, tf_Kac);


x = linspace(-1,1.5);
y = linspace(-0.5,1);
[X,Y] = meshgrid(x,y);
XZ = repmat(X,1,1,3);
YZ = repmat(Y,1,1,3);
ba1 = reshape(A(:,1),1,1,3);
ba2 = reshape(A(:,2),1,1,3);
br = reshape(b,1,1,3);
z = bsxfun(@times,ba1,XZ )+bsxfun(@times,ba2,YZ );
z = bsxfun(@minus,br,z);
z(z<0)=0;
z=0.5*z.^2;
Z=squeeze(sum(z,3));
figure
contour(X,Y,Z,40)
hold on
plot(xA_Kac(1,:),xA_Kac(2,:),'b+')
plot(xA_GS(1,:),xA_GS(2,:),'ro')
plot(xA_In(1,:),xA_In(2,:),'g*')
line([-0.5,1],[-0.5,1]);
line([0,1.5],[1,-0.5]);
line([-1,1.5],[0,0]);

