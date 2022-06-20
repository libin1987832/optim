%�ÿ��������㷨��BBP�����͹���ι滮����
%���  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [x, F]=BBP(A, b, F, maxIt, eps) % epsilon2��ʾ��e��ʾ����ȫΪ1�ĺ�������e0��ʾ��ʾ����ȫΪ0�ĺ�������optionsΪ�˲���� quadprog�м�ļ������
i=0;  %��������
M=[A -e';e 0]; % M����
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
  xk=Z(1:n); %���x
  break
    end
end
w1=find(ZF1<-epsilon2);
v(F)=1;
v(T)=VT;
w2=find(v<epsilon2);
F(w1)=[];%����ָ��F
F=union(F,w2);
i=i+1;  
end   
 end 