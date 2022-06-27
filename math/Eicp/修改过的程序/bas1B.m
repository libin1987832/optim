function [i]=bas1B(A,x1,B,n)
e=ones(n,1);      %����nά��1����
Beta=10^(-5);
Etak=1;
varepsilon=10^(-5);%����
%%  ��ʼѭ��
tic                %��ʼ����ʱ��
k=0;
Etamin=unifrnd (0,1);
Etamax=1/Etamin;
L1=x1'*B*x1;
L3=x1'*A*x1;
f1x1=-(2/L1).*((L3/L1).*B*x1-A*x1); %f1��ʾ����F�ĵ���
dk=zeros(n,1);
while 1
    %% ����dk
    for t=1:n
        if x1(t)<=Beta*f1x1(t)
            dk(t)=-x1(t);
        elseif x1(t)<Etak*f1x1(t)
            dk(t)= -x1(t);
        else
            dk(t)= -Etak*f1x1(t);
        end
    end
    if norm(dk)<=varepsilon   %�жϾ�����ֹ
        break
    end
    %% ����Deltak
    a0=(dk'*A*x1)*L1-(dk'*B*x1)*L3;
    a1=(dk'*A*dk)*L1-(dk'*B*dk)*L3;
    a2=(dk'*A*dk)*(dk'*B*x1)-(dk'*B*dk)*(dk'*A*x1);
    b=sqrt(a1*a1-4*a2*a0);
    Delta1=(-a1-b)/(2*a2);
    Delta2=(-a1+b)/(2*a2); %%���κ�����������
    if Delta2<Delta1
        c=Delta1;
        Delta1= Delta2;
        Delta2=c;
    end
    %%Deltak��ֵ
    F0=((x1+dk)'*A*(x1+dk))/((x1+dk)'*B*(x1+dk));
    F1=((x1+Delta1*dk)'*A*(x1+Delta1*dk))/((x1+Delta1*dk)'*B*(x1+Delta1*dk));
    F2=((x1+Delta2*dk)'*A*(x1+Delta2*dk))/((x1+Delta2*dk)'*B*(x1+Delta2*dk));
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
    %% ����Etak
    sk=x1-x0;
    f1x0=f1x1;
    L1=x1'*B*x1;
    L3=x1'*A*x1;
    f1x1=-(2/L1).*((L3/L1).*B*x1-A*x1); %f1��ʾ����F�ĵ���
    wk=f1x1-f1x0;
    if sk'*wk<=0
        Etak= Etamax;
    else
        Etak= min (Etamax ,max(Etamin,(sk'*sk)/(sk'*wk)));
    end
end
%% ���
toc;                       %��������ʱ��

lamda=(x1'*A*x1)/(x1'*B*x1) %����ֵ
b=A*x1-lamda*B*x1;
h=x1'*b;
norm(dk);
i=k                        %��������
end