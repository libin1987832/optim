function [i]=SSQPI(A,x1,n)
varepsilon=10^(-5);          %精度
varepsilon1=10^(-12);        %子问题的精度
dk=x1;
xk=x1;
i=0;
XZI3=[];
YZI3=[];
%%  开始计算时间  
tic
R=max(abs(eig(A)));    %A的谱半径                                          %用A的1范数来近似A的谱（代价会不会太大了）
e=zeros(5000,1);
e(1)=R+0.01;                                 
%%  开始外循环
while 1        %大循环开始只负责循环不判断
   %% lambda(k)的取值范围
    tao=1-xk'*xk/2;          %tao的值
    v=A*xk;
    qk=(tao+xk'*v)/(2-2*tao);%lambda_k的上界
    lk=v-xk;
    g=find(xk~=0);           %找出x_k中所以有大于零的分量的位置       
    x=xk(g);
    l=lk(g);
    pk=min(l./x);            %lambda_k的下界
    m=zeros(n,1);
    %%  开始内循环，求dk(因为取M为单位矩阵，所以这里不直接调用二次规划函数)
    while  1   %小循环开始只负责循环不判断
     lambda(1)=pk;
      m=lambda(1)*xk-v;
    for t=1:n
        if lambda(1)*xk(t)>=lk(t)
            dk(t)=m(t);
        else
            dk(t)= - xk(t);
        end
    end
    f1=xk'*dk-tao;   
        a=(pk+qk)/2;
        lambda(2)=a;
        m=lambda(2)*xk-v;
        for t=1:n
            if lambda(2)*xk(t)>=lk(t)
                dk(t)=m(t);
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
   YZI3=[YZI3  norm(dk)];
    %% 求 alphak 
    i=i+1;
    XZI3=[XZI3  i];
    e(i)=abs(a);
    c=(dk'*dk)/2;
    sigma=max(e);
    yk=dk+v-lambda*xk;
    z=(2*c+yk'*xk-tao* lambda-tao*sigma)/(dk'*A*dk+2*c*sigma);
    u=2*tao/(tao+sqrt(tao^2+4*c*tao));
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
dlmwrite('XZI3.txt', XZI3, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZI3.txt', YZI3, 'precision', '%5f', 'delimiter', '\t');
toc;
i;
s=toc;
y=A*xk-lambda*xk;
h=y'* xk;
end 