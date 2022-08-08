function [i]=SQPGB(A,x1,B,n,I)
varepsilon=10^(-5);          %精度
xk=x1;
i=0;
XFB2=[];
YFB2=[];
%%  开始计算时间
tic                           
O=inv(B)*A;
R=max(abs(eig(O))); %A的谱半径                                          %用A的1范数来近似A的谱（代价会不会太大了）
e=zeros(5000,1);
e(1)=R+0.01;  
R1=max(abs(eig(A)))+0.01;
%  M=A+R1*I; %A对称非正定
M=A;%A对称正定
%%  开始外循环
while 1        %大循环开始只负责循环不判断
   %% lambda(k)的取值范围
   Bxk=B*xk;
   tao=1-(xk'*Bxk/2);          %tao的值
    v=A*xk;
    Aeq=(Bxk)';
    lb=-xk;
    %%  求dk,直接调用函数 quadprog，(A对称正定时取，M=A.A对称时取M=A+ R*I)
    [dk,fval,exitflag,output,lambda]=quadprog(M,v,[],[],Aeq,tao,lb,[]) ;  
   lambda=(xk'*M*dk+xk'*v+dk'*M*dk+dk'*v)/(xk'*Bxk+dk'*Bxk);
    if  norm(dk)<=varepsilon   %判断大循环终止条件 dk的模
        break
    end
     YFB2=[YFB2   norm(dk)];
    %% 求 alphak 
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
i  ;            %迭代次数
xk;
y=A*xk-lambda.*(Bxk);
h=y'* xk;
norm(dk);
end 