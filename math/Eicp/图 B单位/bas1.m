%%  BAS本程序目的：求出当A为对称（正定）矩阵，B为单位矩阵的特征值互补问题.
function [i]=bas1(A,x1,n)
e=ones(n,1);      %生成n维的1向量
Beta=10^(-5);
Etak=1;
varepsilon=10^(-5);%精度
XZI1=[];
YZI1=[];
%%  开始循环
tic
k=0;
Etamin=unifrnd (0,1);
Etamax=1/Etamin;
L1=x1'*x1;
L3=x1'*A*x1;
f1x1=-(2/L1).*((L3/L1).*x1-A*x1); %f1表示函数F的导数
dk=x1;
while 1
    %% 计算dk
    for t=1:n
        if x1(t)<=Beta*f1x1(t)
            dk(t)=-x1(t);
        elseif x1(t)<Etak*f1x1(t)
            dk(t)= -x1(t);
        else
            dk(t)= -Etak*f1x1(t);
        end
    end
    if norm(dk)<=varepsilon   %判断精度终止
        break
    end
    YZI1=[YZI1   norm(dk)];
    %% 计算Deltak
    a0=(dk'*A*x1)*L1-(dk'*x1)*L3;
    a1=(dk'*A*dk)*L1-(dk'*dk)*L3;
    a2=(dk'*A*dk)*(dk'*x1)-(dk'*dk)*(dk'*A*x1);
    b=sqrt(a1*a1-4*a2*a0);
    Delta1=(-a1-b)/(2*a2);
    Delta2=(-a1+b)/(2*a2); %%二次函数的两个根
    if Delta2<Delta1
        c=Delta1;
        Delta1= Delta2;
        Delta2=c;
    end
    %%Deltak赋值
    F0=((x1+dk)'*A*(x1+dk))/((x1+dk)'*(x1+dk));
    F1=((x1+Delta1*dk)'*A*(x1+Delta1*dk))/((x1+Delta1*dk)'*(x1+Delta1*dk));
    F2=((x1+Delta2*dk)'*A*(x1+Delta2*dk))/((x1+Delta2*dk)'*(x1+Delta2*dk));
    if Delta2>0&&Delta2<=1
        if Delta1>0&&Delta1<=1
            if F2<F0&&F2<F1
                Deltak=Delta2;
            elseif F1<F0&&F1<F2
                Deltak=Delta1;
            else
                Deltak=1;
            end
        elseif F2<F0
            Deltak=Delta2;
        else
            Deltak=1;
        end
    elseif Delta1>0&&Delta1<=1
        if F1<F0
            Deltak=Delta1;
        else
            Deltak=1;
        end
    else
        Deltak=1;
    end
    x0=x1;
    x3=x1+Deltak*dk;
    Muk=1/(e'*x3);
    x1=x3.* Muk;
    k=k+1;
    XZI1=[XZI1   k];
    %% 计算Etak
    sk=x1-x0;
    f1x0=f1x1;
    L1=x1'*x1;
    L3=x1'*A*x1;
    f1x1=-(2/L1).*((L3/L1).*x1-A*x1); %f1表示函数F的导数
    wk=f1x1-f1x0;
    if sk'*wk<=0
        Etak= Etamax;
    else
        Etak= min (Etamax ,max(Etamin,(sk'*sk)/(sk'*wk)));
    end
end
%% 输出
dlmwrite('XZI1.txt', XZI1, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZI1.txt', YZI1, 'precision', '%5f', 'delimiter', '\t');
toc;                       %计算所用时间
i=k;                         %迭代次数
s=toc;
lamda=(x1'*A*x1)/(x1'*x1); %特征值
b=A*x1-lamda*x1;
h=x1'*b;
end