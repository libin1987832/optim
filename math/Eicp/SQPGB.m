function [i]=SQPGB(A,x1,B,n,I)
varepsilon=10^(-5);          %����
xk=x1;
i=0;
XFB2=[];
YFB2=[];
%%  ��ʼ����ʱ��
tic                           
O=inv(B)*A;
R=max(abs(eig(O))); %A���װ뾶                                          %��A��1����������A���ף����ۻ᲻��̫���ˣ�
e=zeros(5000,1);
e(1)=R+0.01;  
R1=max(abs(eig(A)))+0.01;
%  M=A+R1*I; %A�ԳƷ�����
M=A;%A�Գ�����
%%  ��ʼ��ѭ��
while 1        %��ѭ����ʼֻ����ѭ�����ж�
   %% lambda(k)��ȡֵ��Χ
   Bxk=B*xk;
   tao=1-(xk'*Bxk/2);          %tao��ֵ
    v=A*xk;
    Aeq=(Bxk)';
    lb=-xk;
    %%  ��dk,ֱ�ӵ��ú��� quadprog��(A�Գ�����ʱȡ��M=A.A�Գ�ʱȡM=A+ R*I)
    [dk,fval,exitflag,output,lambda]=quadprog(M,v,[],[],Aeq,tao,lb,[]) ;  
   lambda=(xk'*M*dk+xk'*v+dk'*M*dk+dk'*v)/(xk'*Bxk+dk'*Bxk);
    if  norm(dk)<=varepsilon   %�жϴ�ѭ����ֹ���� dk��ģ
        break
    end
     YFB2=[YFB2   norm(dk)];
    %% �� alphak 
    i=i+1;
    XFB2=[XFB2   i];
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
dlmwrite('XZB2.txt', XFB2, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZB2.txt', YFB2, 'precision', '%5f', 'delimiter', '\t');
toc;
i  ;            %��������
xk;
y=A*xk-lambda.*(Bxk);
h=y'* xk;
norm(dk);
end 