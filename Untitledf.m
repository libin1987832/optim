[f,x]=fminunc('sqrt(((x(1)-x(4))^2+(x(2)+x(3))^2)*((x(1)+x(4))^2+(x(2)-x(3))^2))-x(1)^2-x(3)^2+x(2)^2+x(4)^2',[-1;-10;23;33])
A = rand(3);
x = linspace(-1, 1,10); y = x; z = x;
[XX,YY,ZZ] = meshgrid(x,y,z);                        %// step 1
XYZ = [XX(:) YY(:) ZZ(:)];                           %// step 2
XYZ2 = bsxfun(@times, XYZ, permute(XYZ, [1 3 2]));   %// step 3
f = reshape(XYZ2,[],numel(A))*A(:);                  %// step 4
f = reshape(f, size(XX));                            %// step 5
level = 1;
s = isosurface(XX,YY,ZZ,f,level);
patch(s, 'EdgeColor','none','FaceColor','blue');     %// st
% for i =1 :1000
% A=rand(4,4)*10+123;
% ATA = A'*A;
% l1 = max(eig(ATA));
% B = rand(1,4);
% B = round(B);
% l2 = max(eig(ATA*diag(B)));
% if l1 > l2
%     save('ss','A','B')
%     [l1 l2];
% end
% end
% 
% A =  [3 1;1 6];
% 
% C=A'*A*[1 0;0 0];
% [max(eig(C)) max(eig(A'*A))]

+%% plot w b in one dimension
+% x = linspace(-0.1,0.1);
+% y = linspace(-0.1,0.3);
+% [X,Y] = meshgrid(x,y);
+% XY = permute([X(:) Y(:)], [3 2 1]);
+% % elementary wise mutilp
+% XYsum = bsxfun(@times,[A1;B1], XY);
+% XYs = sum(XYsum,2);
+% XYnorm = reshape(XYs,[],numel(X));
+% XYnorm = ones(size(XYnorm))-XYnorm;
+% XYnorm(XYnorm<0) = 0;
+% XYn = sum(XYnorm.*XYnorm,1);
+% Z = reshape(XYn,size(X));
+% [m,nn] = size(X);
+% [rows,columns] = find(Z==min(Z(:))); 
+% minbn =10;