% 罚函数方法
% （海森矩阵+防止错误项）* 方向 = 梯度
%  搜索步长
function [x,fk]=GNP(x0,M,delt,e,A,b)
    k=0;  
    d1=det1F(x0,A,b,M);
    while(norm(d1,inf)>e)
        AA=det2F(x0,A,b,M)+diag(ones(size(A',1),1)).*delt;
        d1=det1F(x0,A,b,M);
        % 计算方向
        p=-1*AA\d1;
        % 计算步长
        a=fsearchaM(A,b,x0,p,M);
        x0=x0+a*p;
        k=k+1;
        fQ(A,b,x0,M)
    end
    x0(x0<0)=0;
    x=x0;
    fk=fq(A,b,x);