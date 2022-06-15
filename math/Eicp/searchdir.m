function d = searchdir(A, B, M, x, method)
    if method == 1
    Bxk=B*xk;
    tao=1-(xk'*Bxk/2);          %tao的值
    v=A*xk;
    Aeq=(Bxk)';
    lb=-x;
    %%  求dk,直接调用函数 quadprog，(A对称正定时取，M=A.A对称时取M=A+ R*I)
    [dk,fval,exitflag,output,lambda]=quadprog(M,v,[],[],Aeq,tao,lb,[]) ;  