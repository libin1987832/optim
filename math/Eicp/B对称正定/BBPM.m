%�ÿ��������㷨��BBP�����͹���ι滮����
%���  min f(x)=1/2x'*A*x+yk'*x
%      st. e'*x=1
%          x>=0
function [xk]=BBPM(A,Bxk,tao,v,n,xk1,epsilon2) % epsilon2��ʾ��e��ʾ����ȫΪ1�ĺ�������e0��ʾ��ʾ����ȫΪ0�ĺ�������optionsΪ�˲���� quadprog�м�ļ������
i=0;  %��������
M=[A -Bxk;Bxk' 0]; % M����
% Z=zeros(n+1,1);
h=[v;-tao];
Z=h;
N=1:(n+1);
% N=N';
F=N;
v1=F;
xk1=[xk1;0];
while 1
 T=setdiff(N, F);
 MFF=M(F,F);
 MTF=M(T,F);
 hF=h(F);
 hT=h(T);
ZF=MFF\(-hF);
VT=hT+MTF*ZF;
xkF=xk1(F);
xkF1=xkF(1:end-1);
ZF1=ZF(1:end-1);
if ZF1>=-xkF1
    if VT>=-epsilon2
  Z(F)=ZF;
  Z(T)=-xkF(T);    
  xk=Z(1:n); %���x
  break
    end
end
w1=find(ZF1<-epsilon2);
v1(F)=1;
v1(T)=VT;
w2=find(v1<epsilon2);
F(w1)=[];%����ָ��F
F=union(F,w2);
i=i+1;  
end   
 end 