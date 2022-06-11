%%  SSQP��I��������Ŀ�ģ������AΪ�Գƣ�����������BΪ�Գ��������������ֵ��������.M=I
function [i]=SSQPIB(A,x1,B,n)
varepsilon=10^(-5);          %����
varepsilon1=10^(-12);        %������ľ���
dk=x1;
xk=x1;
i=0;
XFB3=[];
YFB3=[];
%%  ��ʼ����ʱ��
tic      
O=inv(B)*A;
R=max(abs(eig(O)));    %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;                                 
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   %% lambda(k)��ȡֵ��Χ
    Bxk=B*xk;
    tao=1-(xk'*Bxk/2);          %tao��ֵ
    v=A*xk;
    qk=(tao+Bxk'*v)/(Bxk'*Bxk);%lambda_k���Ͻ�
    lk=v-xk;
    g=find(Bxk>0);  %�ҳ�Bx_k�������д�����ķ�����λ��       
    x=Bxk(g);
    l=lk(g);
    pk=min(l./x);            %lambda_k���½�
    %%  ��ʼ��ѭ������dk(��ΪȡMΪ��λ�����������ﲻֱ�ӵ��ö��ι滮����)
    while  1   %Сѭ����ʼֻ����ѭ�����ж�
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
        if qk-pk<varepsilon1  %�ж�Сѭ����ֹ����
            break 
        end
    end
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
     YFB3=[YFB3   norm(dk)];
    %% �� alphak 
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
i   ;           %��������
lambda=a;
y=A*xk-lambda.*(Bxk);
h=y'* xk;
norm(dk);
end 