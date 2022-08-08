%%  SSQP(D)本程序目的：求出当A为对称（正定）矩阵，B为单位矩阵的特征值互补问题. M=diag(A+R*I)/M=diag(A)
function [i]=SSQPD(A,x1,n)
varepsilon=10^(-5);          %精度
varepsilon1=10^(-12);        %子问题的精度
dk=x1;
xk=x1;
i=0;
%%  开始计算时间  
R=max(abs(eig(A)));    %A的谱半径                                          %用A的1范数来近似A的谱（代价会不会太大了）
e=zeros(5000,1);
e(1)=R+0.01;
I=eye(n);
% M=A;%A对称正定
M=A+(R+0.01)*I;%A对称非正定
Q=diag(M);
Q1=1./Q;
Q2=diag(Q1);
tic 
%%  开始外循环
while 1        %大循环开始只负责循环不判断
   %% lambda(k)的取值范围
    tao=1-xk'*xk/2;          %tao的值
    v=A*xk;
    qk=(tao+xk'*Q2*v)/(xk'*Q2*xk);%lambda_k的上界
    lk=v-Q.*xk;
    g=find(xk~=0);           %找出x_k中所以有大于零的分量的位置       
    x=xk(g);
    l=lk(g);
    pk=min(l./x);            %lambda_k的下界
    m=zeros(n,1);
    %%  开始内循环，求dk
    while  1   %小循环开始只负责循环不判断
     lambda(1)=pk;
      m=lambda(1)*xk-v;
    for t=1:n
        if lambda(1)*xk(t)>lk(t)
            dk(t)=m(t)*Q1(t);
        else
            dk(t)= - xk(t);
        end
    end
    f1=xk'*dk-tao;   
        a=(pk+qk)/2;
        lambda(2)=a;
        m=lambda(2)*xk-v;
        for t=1:n
            if lambda(2)*xk(t)>lk(t)
                dk(t)=m(t)*Q1(t);
            else
                dk(t)=- xk(t);
            end
        end
        f2=xk'*dk-tao;
        if  f1*f2<=0
            qk=a;
        else
            pk=a;
        end
        if qk-pk<varepsilon1  %判断小循环终止条件
            break 
        end
    end
    lambda=a;
    if  norm(dk)<=varepsilon   %判断大循环终止条件 dk的模
        break
    end
    %% 求 alphak 
    i=i+1;
    e(i)=abs(a)+0.1;
    c=(dk'*dk)/2;
    sigma=max(e);
    z=-(v'*dk+tao*sigma)/(dk'*A*dk+2*c*sigma);
    u=(-tao+sqrt(tao^2+4*c*tao))/(2*c);
   if tao>0
            if z>=1
               alphak=1;
            elseif z<1&&z>u
                    alphak=z;
            else
                alphak=u;
            end
            
   else
            if z>=1
                alphak=1;
            else
                alphak=z;
            end
    end
    xk=xk+alphak*dk;
end
toc;
xk;
lambda
i

y=A*xk-lambda*xk;
h=y'* xk;
end 