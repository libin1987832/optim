% ��ȷ������⣨�м�ʹ��MATLAB�Դ����㷨����Լ������С�������⣩
% ���Ǹ���С��������ת��� ��ͨ����С��������?
function [x0,f0]=alg2(A,b,x0,e)
[m,n]=size(x0);
options = optimoptions('LSQLIN'); 
 options.OptimalityTolerance=e;
index=0;
y0=b-A*x0;
y0(y0<0)=0;
res0=norm(min(x0,A'*A*y0));
t=tic;
while 1
    y0=b-A*x0;
    y0(y0<0)=0;
    %     z0=A*x0-b;
    %     z0(z0<0)=0;
    bk=y0+A*x0;
	 %����ʹ��MATLAB�Դ����㷨
%     [x1,f1,residual,exitflag]=lsqlin(A,bk,[],[],[],[],zeros(m,1),Inf*ones(m,1),x0,options);
   [x1,f1,residual,exitflag]=lsqlin(A,bk,[],[],[],[],[],[],x0,options);
    
    f0=b-A*x0;
    f0(f0<0)=0;
%     f0=0.5*norm(f0);
    f0=0.5*f0'*f0;
    x0=x1;
    
    res1=norm(min(x0,A'*A*y0));
     fprintf('index:%d,exit %d,f:%f,res1:%f,res0:%f,ratio:%f!\n',index,exitflag,f1,res1,res0,res1/res0);
    index=index+1;
    if res1/res0<e
        break;
    end
%     if norm(f1-f0)<e
%         break;
%     end
end
toc(t);
% fprintf('exact run time cost:%f\n',etime(clock,t));
