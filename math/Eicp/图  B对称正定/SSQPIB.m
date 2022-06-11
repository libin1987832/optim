%%  SSQP（I）本程序目的：求出当A为对称（正定）矩阵，B为对称正定矩阵的特征值互补问题.M=I
function [i]=SSQPIB(A,x1,B,n)
varepsilon=10^(-5);          %精度
varepsilon1=10^(-12);        %子问题的精度
dk=x1;
xk=x1;
i=0;
XFB3=[];
YFB3=[];
%%  开始计算时间
tic      
O=inv(B)*A;
R=max(abs(eig(O)));    %A的谱半径                                          %用A的1范数来近似A的谱（代价会不会太大了）
e=zeros(5000,1);
e(1)=R+0.01;                                 
%%  开始外循环
while 1        %大循环开始只负责循环不判断
   %% lambda(k)的取值范围
    Bxk=B*xk;
    tao=1-(xk'*Bxk/2);          %tao的值
    v=A*xk;
    qk=(tao+Bxk'*v)/(Bxk'*Bxk);%lambda_k的上界
    lk=v-xk;
    g=find(Bxk>0);  %找出Bx_k中所以有大于零的分量的位置       
    x=Bxk(g);
    l=lk(g);
    pk=min(l./x);            %lambda_k的下界
    %%  开始内循环，求dk(因为取M为单位矩阵，所以这里不直接调用二次规划函数)
    while  1   %小循环开始只负责循环不判断
     lambda1=pk;
      m=lambda1*Bxk-v;
    for t=1:n
        if lambda1*Bxk(t)>lk(t)
            dk(t)=m(t);
        else
            dk(t)= - xk(t);
        end
    end
    f1=Bxk'*dk-tao;   
        a=(pk+qk)/2;
        lambda2=a;
        m=lambda2*Bxk-v;
        for t=1:n
            if lambda2*Bxk(t)>lk(t)
                dk(t)=m(t);
            else
                dk(t)=- xk(t);
            end
        end
        f2=Bxk'*dk-tao;
        if  f1*f2<=0
            qk=a;
        else
            pk=a;
        end
        if qk-pk<varepsilon1  %判断小循环终止条件
            break 
        end
    end
    if  norm(dk)<=varepsilon   %判断大循环终止条件 dk的模
        break
    end
     YFB3=[YFB3   norm(dk)];
    %% 求 alphak 
    i=i+1;
    XFB3=[XFB3  i];
    e(i)=abs(a)+0.1;
    c=(dk'*B*dk)/2;
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
dlmwrite('XZB3.txt', XFB3, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZB3.txt', YFB3, 'precision', '%5f', 'delimiter', '\t');
toc;
i   ;           %迭代次数
lambda=a;
y=A*xk-lambda.*(Bxk);
h=y'* xk;
norm(dk);
end 