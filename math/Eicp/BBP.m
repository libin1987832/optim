%用块主枢轴算法（BBP）求解凸二次规划问题
%求解  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [x, F]=BBP(A, b, F, maxIt, eps) % epsilon2表示误差，e表示分量全为1的横向量，e0表示表示分量全为0的横向量，options为了不输出 quadprog中间的计算过程
i=0;  %迭代次数
M=[A -e';e 0]; % M矩阵
Z=zeros(n+1,1);
h=[yk;-1];
N=1:(n+1);
% N=N';
F=N;
v=F;
while 1
 T=setdiff(N, F);
 MFF=M(F,F);
 MTF=M(T,F);
 hF=h(F);
 hT=h(T);
ZF=MFF\(-hF);
VT=hT+MTF*ZF;
ZF1=ZF(1:end-1);
if ZF1>=-epsilon2
    if VT>=-epsilon2
  Z(F)=ZF;
  Z(T)=0;    
  xk=Z(1:n); %输出x
  break
    end
end
w1=find(ZF1<-epsilon2);
v(F)=1;
v(T)=VT;
w2=find(v<epsilon2);
F(w1)=[];%更新指标F
F=union(F,w2);
i=i+1;  
end   
 end 