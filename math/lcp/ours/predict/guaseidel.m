function [xk,n,xA] = guaseidel(A,b,x0,eps,it_max)
% 求线性方程组的Gauss-Seidel迭代法，调用格式为

%  [x, k] = guaseidel(A,b,x0,eps,it_max)

%  其中, A 为线性方程组的系数矩阵，b 为常数项，eps 为精度要求，默认为1e-5,

%  it_max 为最大迭代次数，默认为100

%  x 为线性方程组的解，k迭代次数
if nargin == 3
    eps = 1.0e-6
    it_max= 200
elseif nargin == 4
    it_amx = 200
elseif nargin <3
    disp('输入参数个数不足3个');
    return;
end
xA=[];
% D = diag(diag(A));%求A的对角矩阵
% L = -tril(A,-1);%求A的下三角矩阵,不带对角线
% U = -triu(A,1);%求A的上三角矩阵
% G = (D-L)\U;
% f = (D-L)\b;
% x = G*x0+f;
[ms,ns]=size(A);
n=1;  %迭代次数
s=1;
%while norm(xk-x0)>=eps
%     x0 = x;
%     x = G*x0+f;
for i=1:it_max
    xk=zeros(ns,1);
    for j=1:ns
        sd=s*(b(j)-A(j,:)*[xk(1:(j-1));x0(j:ns)])/A(j,j);
        xk(j)=x0(j)-sd;
    end
    x0=xk;
    xA=[xA xk];
    n = n+1;
    if(n>=it_max)
        disp('Warning:迭代次数太多,可能不收敛');
        return;
    end
end

