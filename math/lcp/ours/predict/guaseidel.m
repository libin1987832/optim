function [x,n,xA] = guaseidel(A,b,x0,eps,it_max)
% �����Է������Gauss-Seidel�����������ø�ʽΪ
 
%  [x, k] = guaseidel(A,b,x0,eps,it_max)
 
%  ����, A Ϊ���Է������ϵ������b Ϊ�����eps Ϊ����Ҫ��Ĭ��Ϊ1e-5,
 
%  it_max Ϊ������������Ĭ��Ϊ100
 
%  x Ϊ���Է�����Ľ⣬k��������
if nargin == 3
    eps = 1.0e-6
    it_max= 200
elseif nargin == 4
    it_amx = 200
elseif nargin <3
    disp('���������������3��');
    return;
end
xA=[];
D = diag(diag(A));%��A�ĶԽǾ���
L = -tril(A,-1);%��A�������Ǿ���,�����Խ���
U = -triu(A,1);%��A�������Ǿ���
G = (D-L)\U;
f = (D-L)\b;
x = G*x0+f;
n=1;  %��������
while norm(x-x0)>=eps
    x0 = x;
    x = G*x0+f;
    xA=[xA x];
    n = n+1;
    if(n>=it_max)
        disp('Warning:��������̫��,���ܲ�����');
        return;
    end
end
 
