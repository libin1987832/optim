function [i]=SQPGB(A,x1,B,n,I)
varepsilon=10^(-5);          %����
options = optimoptions('quadprog','Display','off');%Ϊ�˲���� quadprog�м�Ĺ���
xk=x1;
i=0;                          
O=inv(B)*A;
R=max(abs(eig(O))); %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;  
R1=max(abs(eig(A)))+0.01;

M=A+R1*I; %A�ԳƷ�����
%  M=A;%A�Գ�����
 M=(M+M')/2;
%%  ��ʼ����ʱ��
 tic 
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   %% lambda(k)��ȡֵ��Χ
   Bxk=B*xk;
   tao=1-(xk'*Bxk/2);          %tao��ֵ
    v=A*xk;
    Aeq=(Bxk)';
    lb=-xk;
    %%  ��dk,ֱ�ӵ��ú��� quadprog��(A�Գ�����ʱȡ��M=A.A�Գ�ʱȡM=A+ R*I)
    
    [dk,~,~,~,~]=quadprog(M,v,[],[],Aeq,tao,lb,[],[],options) ;  
   lambda=(xk'*M*dk+xk'*v+dk'*M*dk+dk'*v)/(xk'*Bxk+dk'*Bxk);
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
    %% �� alphak 
    i=i+1;
    e(i)=abs(lambda)+0.01;
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
toc;

lambda

i  %��������

end 