%%  SSQP��D��������Ŀ�ģ������AΪ�Գƣ�����������BΪ�Գ��������������ֵ��������.M=diag(A+R*B)/M=A
function [i]=SSQPDB(A,x1,B,n,I)
varepsilon=10^(-5);          %����
varepsilon1=10^(-12);        %������ľ���
dk=x1;
xk=x1;
i=1;
O=inv(B)*A;
R=max(abs(eig(O)));    %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;
R1=max(abs(eig(A)))+0.01;
%  M=A+R1*I;%A�ԳƷ�����
 M=A; %A�Գ�����
Q=diag(M);
Q1=1./Q;
Q2=diag(Q1);
%%  ��ʼ����ʱ��
tic 
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   %% lambda(k)��ȡֵ��Χ
    Bxk=B*xk;
   tao=1-xk'*Bxk/2;          %tao��ֵ
    v=A*xk;
    qk=(tao+Bxk'*Q2*v)/(Bxk'*Q2*Bxk);%lambda_k���Ͻ�
    lk=v-Q.*xk;
    g=find(Bxk>0);           %�ҳ�Bx_k�������д�����ķ�����λ��       
    x=Bxk(g);
    l=lk(g);
    pk=min(l./x) ;          %lambda_k���½�
    %%  ��ʼ��ѭ������dk(��ΪȡMΪ��λ�����������ﲻֱ�ӵ��ö��ι滮����)
    while  1   %Сѭ����ʼֻ����ѭ�����ж�
     lambda(1)=pk;
      m=lambda(1)*Bxk-v;
     
      
    for t=1:n
        if lambda(1)*Bxk(t)>lk(t)
            dk(t)=m(t)*Q1(t);
        else
            dk(t)= - xk(t);
        end
    end
    f1=Bxk'*dk-tao;   
        a=(pk+qk)/2;
        lambda(2)=a;
        m=lambda(2)*Bxk-v;
        
        for t=1:n
            if lambda(2)*Bxk(t)>lk(t)
                dk(t)=m(t)*Q1(t);
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
    lambda=a;
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
    %% �� alphak 
    i=i+1;
    e(i)=abs(a)+0.1;
    c=(dk'*B*dk)/2;
    sigma=max(e);
    z=-(v'*dk+tao*sigma)/(dk'*A*dk+2*c*sigma);
    u=(-tao+sqrt(tao^2+4*c*tao))/(2*c);
   if tao>0
            if z>=1
               alphak=1;
            elseif z>u
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

i  %��������
y=A*xk-lambda.*(Bxk);
h=y'* xk;
norm(dk);

end