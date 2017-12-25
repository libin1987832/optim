function [x0,f0]=alg1(A,b,x0,e)
[m,n]=size(x0);
options = optimoptions('LSQLIN'); 
options.OptimalityTolerance=e;
index=0;
while 1
    y0=b-A*x0;
    y0(y0<0)=0;
    %     z0=A*x0-b;
    %     z0(z0<0)=0;
    bk=y0+A*x0; 
    [x1,f1,residual,exitflag]=lsqlin(A,bk,[],[],[],[],zeros(m,1),Inf*ones(m,1),x0,options);
    fprintf('index:%d,exit %d!,%f\n',index,exitflag,residual);
    index=index+1;
    f0=b-A*x0;
    f0(f0<0)=0;
    f0=0.5*norm(f0);
    x0=x1;
    if norm(f1-f0)<e
        break;
    end
end