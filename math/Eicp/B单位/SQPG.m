function [i]=SQPG(A,x1,n)
varepsilon=10^(-5);   %����
options = optimoptions('quadprog','Display','off');%Ϊ�˲���� quadprog�м�Ĺ���
xk=x1;
i=0;
%%  ��ʼ����ʱ��                          
R=max(abs(eig(A))); %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;    
I=eye(n);
% M=A;       %A�Գ�����
M=A+(R+0.01)*I ;  %A�ԳƷ�����
M=(M+M')/2;
tic 
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   tao=1-(xk'*xk/2);          %tao��ֵ
    v=A*xk;
    Aeq=(xk)';
    lb=-xk;
    %%  ��dk,ֱ�ӵ��ú��� quadprog��(A�Գ�����ʱȡ��M=A.A�Գ�ʱȡM=A+ R*I)
%     [dk,~,~,~,~]=quadprog(M,v,[],[],Aeq,tao,lb,[],[],options) ;
        dk=quadprog(M,v,[],[],Aeq,tao,lb,[],[],options) ;
   lambda=(xk'*M*dk+xk'*v+dk'*M*dk+dk'*v)/(xk'*xk+dk'*xk);
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
    %% ��alphak 
    i=i+1;
   e(i)=abs(lambda);
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
i %��������
s=toc;%����ʱ��
y=A*xk-lambda*xk;
h=y'* xk;

end