function x=GNP(x0,M,delt,e,A,b)
%������
    k=0;  
    d1=det1F(x0,A,b,M);
    while(d1>e)
        AA=det2F(x0,A,b,M)+diag(size(A',1)).*delt;
        d1=det1F(x0,A,b,M);
        p=-1*d1/AA;
        a=fsearchaM(A,b,x0,p,M)
        x0=x0+a*p;
        k=k+1;
    end
    x0(x0<0)=0;
    x=x0;