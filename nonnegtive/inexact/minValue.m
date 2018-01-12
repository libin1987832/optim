function [minNum,v]=minValue(x,scale,number,A,b,f,c)
[m,n]=size(x);
fx=f(A,b,x);
v=zeros(m,number);
minNum=0;
index=0;
while 1
    %随机产生数字
    xi=x+(2*rand(m,n)-ones(m,n))*scale;
    g=c(xi);
    %满足约束条件
    if g>0
        fxi=f(A,b,xi);
        %比给定最小值还要小
        if fxi<fx
            minNum=minNum+1;
            v(:,minNum)=xi;
        end
        %只要符合约束的值才能算有效值
        index=index+1;
        %实验值比给定的要多 就停止判断
        if index>number
            break;
        end
    end
end; 