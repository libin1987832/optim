function d = searchdir(A, B, M, x, method)
    if method == 1
    Bxk=B*xk;
    tao=1-(xk'*Bxk/2);          %tao��ֵ
    v=A*xk;
    Aeq=(Bxk)';
    lb=-x;
    %%  ��dk,ֱ�ӵ��ú��� quadprog��(A�Գ�����ʱȡ��M=A.A�Գ�ʱȡM=A+ R*I)
    [dk,fval,exitflag,output,lambda]=quadprog(M,v,[],[],Aeq,tao,lb,[]) ;  