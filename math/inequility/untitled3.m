A=[1 2;3 4];
b=[3;4];
xs = [-2;2.5];
x=[2;0];
ax = b-A*x;
ax(ax<0) =0 ;
I=eye(2);
-(x-xs)'*(I+A'*A)*A'*ax
d=A'*ax;
dd = d'*d