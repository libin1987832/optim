% 罚函数方法
% （海森矩阵+防止错误项）* 方向 = 梯度
%  搜索步长
function [x,fk]=GNP(x0,M,delt,e,A,b)
    k=0;  
    d1=det1F(x0,A,b,M);
    y0=b-A*x0;
    y0(y0<=0)=0;
    res0=norm(min(x0,A'*A*y0));
    t=tic;
%     while(norm(d1,inf)>e)
    while(1)
        AA=det2F(x0,A,b,M)+diag(ones(size(A',1),1)).*delt;
        d1=det1F(x0,A,b,M);
        % 计算方向
        p=-1*AA\d1;
        % 计算步长
        a=fsearchaM(A,b,x0,p,M);
        x0=x0+a*p;
        
        y0=b-A*x0;
        y0(y0<=0)=0;
        res1=norm(min(x0,A'*A*y0));
       
        f=fQ(A,b,x0,M);
        fprintf('index:%d,f1:%f,res1:%f,res0:%f,ratio:%f\n',k,f,res1,res0,res1/res0);
         k=k+1;
        if res1/res0<e
            break;
        end
    end
    toc(t);
    x0(x0<0)=0;
    x=x0;
    fk=fq(A,b,x);