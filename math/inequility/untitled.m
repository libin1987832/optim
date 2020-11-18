clear;clc

A=[1;2];
B=[-1;0;4];
A1=[A [-1;-1]];
B1=[-B [1;1;1]];
b1=ones(5,1);
maxIter = 20;
x0=zeros(2,1);
[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,[A1;B1],b1,maxIter);
xkh1


x = linspace(0,0.5);
y = linspace(0,0.5);
[X,Y] = meshgrid(x,y);
XY = permute([X(:) Y(:)], [3 2 1]);
XYsum = bsxfun(@times,[A1;B1], XY);
XYs = sum(XYsum,2);
XYnorm=reshape(XYs,[],numel(X));
XYnorm(XYnorm<0)=0;
Z=zeros(size(X));

  %  contour(X,Y,Z)

% A = rand(3);
% x = linspace(-1, 1,2); y = x; z = x;
% [XX,YY,ZZ] = meshgrid(x,y,z);                        %// step 1
% XYZ = [XX(:) YY(:) ZZ(:)];                           %// step 2
% XYZ2 = bsxfun(@times, XYZ, permute(XYZ, [1 3 2]));   %// step 3
% f = reshape(XYZ2,[],numel(A))*A(:);                  %// step 4
% f = reshape(f, size(XX));                            %// step 5
% level = 1;
% s = isosurface(XX,YY,ZZ,f,level);
% patch(s, 'EdgeColor','none','FaceColor','blue');  

% gamm = 0;
% px1 = 0:0.2:1;
% p1 = [px1' -px1'+1.5];fm1 = size(p1,1);
% p2 = [px1' -px1'+0.5];fm2 = size(p2,1);
% %plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');
% A1(1:fm1,:) = p1;
% A1(fm1+1:fm1+fm2,:) = -p2;
% b1 = [(1 + gamm)*ones(fm1,1);(1-gamm)*ones(fm2,1)];% 
% x0=zeros(size(A1,2),1);maxIter = 10;
% [xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
% ww = [ xkh1';xkh1'; ];bb = [ gamm+1 ;gamm-1;];
% d = lineData(ww , bb, [0,1], [0,1]);% »­¶àÌõÏß
% %hold on 
% %line(d(:,[1,2])',d(:,[3,4])')
% [[xkh1' gamm+1]/norm(xkh1) norm(rkh);[xkh1' 1-gamm]/norm(xkh1) norm(A1'*rkh)]
