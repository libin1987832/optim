[f,x]=fminunc('sqrt(((x(1)-x(4))^2+(x(2)+x(3))^2)*((x(1)+x(4))^2+(x(2)-x(3))^2))-x(1)^2-x(3)^2+x(2)^2+x(4)^2',[-1;-10;23;33])
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