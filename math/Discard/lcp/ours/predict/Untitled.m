C=[4,-3;-3,8];
q=[5;-16];
% -[0,3/4;0,9/32]^2*[5/4;-49/32]-[0,3/4;0,9/32]*[5/4;-49/32]-[5/4;-49/32];
 R=[0,3/4;0,9/32];
 b=-[5/4;-49/32];
  D=tril(C);
 DR=diag(ones(2,1))-R;
 DRi=inv(DR);
 R^4*[1;1]+R^3*b+R^2*b+R*b+b;
 R^3*b+R^2*b+R*b+b
 DRi*(diag(ones(2,1))-R^4)*b
 DRi*b
 C\q
%  (D*b)\C
%  q\C
%  x0=[1;1];
%    [xks,ress]=splitS(C,q,1,x0,1)
%    b+R*x0
%  b
% R*b
% R*b+b
% R*(R*b+b)
% R*(R*b+b)+b
% R*(R*(R*b+b)+b)
% R*(R*(R*b+b)+b)+b