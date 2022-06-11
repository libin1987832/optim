%%  SSQP(D)������Ŀ�ģ������AΪ�Գƣ�����������BΪ��λ���������ֵ��������. M=diag(A+R*I)/M=diag(A)
function [i]=SSQPD(A,x1,n)
varepsilon=10^(-5);          %����
varepsilon1=10^(-12);        %������ľ���
dk=x1;
xk=x1;
i=0;
%%  ��ʼ����ʱ��  
R=max(abs(eig(A)));    %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;
I=eye(n);
% M=A;%A�Գ�����
M=A+(R+0.01)*I;%A�ԳƷ�����
Q=diag(M);
Q1=1./Q;
Q2=diag(Q1);
tic 
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   %% lambda(k)��ȡֵ��Χ
    tao=1-xk'*xk/2;          %tao��ֵ
    v=A*xk;
    qk=(tao+xk'*Q2*v)/(xk'*Q2*xk);%lambda_k���Ͻ�
    lk=v-Q.*xk;
    g=find(xk~=0);           %�ҳ�x_k�������д�����ķ�����λ��       
    x=xk(g);
    l=lk(g);
    pk=min(l./x);            %lambda_k���½�
    m=zeros(n,1);
    %%  ��ʼ��ѭ������dk
    while  1   %Сѭ����ʼֻ����ѭ�����ж�
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
        if qk-pk<varepsilon1  %�ж�Сѭ����ֹ����
            break 
        end
    end
    lambda=a;
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
    %% �� alphak 
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