% 罚函数的根据方向求步长搜索算法
% 满足wolf准则即表示步长可以找到
function  a=fsearchaM(A,b,x0,p,M)
% 参数设置
t=0.5;
a0=1;
c1=0.01;
c2=0.9;
q1=fQ(A,b,x0+a0*p,M);
q0=fQ(A,b,x0,M);

g1=det1F(x0+a0*p,A,b,M);
g0=det1F(x0,A,b,M);
% wolf准则
%while x0+a0*p<0||q1>=q0+c1*a0*g0'*p||g1'*p<=c2*g0'*p;
% while q1>=q0+c1*a0*g0'*p||g1'*p<=c2*g0'*p;
while q1>=q0+c1*a0*g0'*p   
    % 减少步长，尝试是否满足wolf条件
    t=t/2;
    a0=a0*t;
    % 更新新步长下函数值
    q1=fQ(A,b,x0+a0*p,M);
    % 更新梯度
    g0=det1F(x0+a0*p,A,b,M);
%     q1=fq(A,b,x0+a0*p);
%  
%     g1=fdetq(A,b,x0+a0*p);
end
a=a0;
end

