


%---------
% A=[1 2;3 4];
% b=[3;4];
% xs = [-2;2.5];
% x=[2;0];
% axt = b-A*x;
% ax=axt;
% ax(ax<0) =0 ;
% I=eye(2);
% d=A'*ax;
% -(x-xs)'*(I+A'*A)*d
% dd = d'*d
% 
% z = -axt;
% z(z<0)=0;
% z'*A*d
% 
% dx=A*(x-xs)-(z);
% dx'*dx