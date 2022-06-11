function [i]=SSQPI(A,x1,n)
varepsilon=10^(-5);          %����
varepsilon1=10^(-12);        %������ľ���
dk=x1;
xk=x1;
i=0;
XZI3=[];
YZI3=[];
%%  ��ʼ����ʱ��  
tic
R=max(abs(eig(A)));    %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;                                 
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   %% lambda(k)��ȡֵ��Χ
    tao=1-xk'*xk/2;          %tao��ֵ
    v=A*xk;
    qk=(tao+xk'*v)/(2-2*tao);%lambda_k���Ͻ�
    lk=v-xk;
    g=find(xk~=0);           %�ҳ�x_k�������д�����ķ�����λ��       
    x=xk(g);
    l=lk(g);
    pk=min(l./x);            %lambda_k���½�
    m=zeros(n,1);
    %%  ��ʼ��ѭ������dk(��ΪȡMΪ��λ�����������ﲻֱ�ӵ��ö��ι滮����)
    while  1   %Сѭ����ʼֻ����ѭ�����ж�
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
        if qk-pk<varepsilon1  %�ж�Сѭ����ֹ����
            break 
        end
    end
    lambda=a;
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
   YZI3=[YZI3  norm(dk)];
    %% �� alphak 
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