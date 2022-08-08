function [i]=SQPG(A,x1,n)
varepsilon=10^(-5);   %精度
xk=x1;
dk=x1;
i=0;
XZI2=[];
YZI2=[];
%%  开始计算时间
tic                           
R=max(abs(eig(A))); %A的谱半径                                          %用A的1范数来近似A的谱（代价会不会太大了）
e=zeros(5000,1);
e(1)=R+0.01;    
I=eye(n);
 M=A;       %A对称正定
% M=A+(R+0.01)*I ;  %A对称非正定
M=(M+M')/2;
%%  开始外循环
while 1        %大循环开始只负责循环不判断
   tao=1-(xk'*xk/2);          %tao的值
    v=A*xk;
    Aeq=(xk)';
    lb=-xk;
    %%  求dk,直接调用函数 quadprog，(A对称正定时取，M=A.A对称时取M=A+ R*I)
    [dk,fval,exitflag,output,lambda]=quadprog(M,v,[],[],Aeq,tao,lb,[]) ;  
   lambda=(xk'*M*dk+xk'*v+dk'*M*dk+dk'*v)/(xk'*xk+dk'*xk);
    if  norm(dk)<=varepsilon   %判断大循环终止条件 dk的模
        break
    end
      YZI2=[YZI2   norm(dk)];
    %% 求alphak 
    i=i+1;
    XZI2=[XZI2  i];
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
dlmwrite('XZI2.txt', XZI2, 'precision', '%5f', 'delimiter', '\t');
dlmwrite('YZI2.txt', YZI2, 'precision', '%5f', 'delimiter', '\t');
toc;
i;%迭代次数
s=toc;%迭代时间
y=A*xk-lambda*xk;
h=y'* xk;
end