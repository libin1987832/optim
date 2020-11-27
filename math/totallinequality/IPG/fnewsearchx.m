%IPG ¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿? ¿¿¿¿¿¿¿¿¿¿
function x=fnewsearchx(A,b,x0,e,det)
%Ê¹ÓÃ²»ËÑË÷²½³¤·½·¨£¬ÕÒ×îÓÅx
g=fdetq(A,b,x0);
i=0;
%while norm(x0.*g,inf)>e||min(g)<-e
while norm(x0.*g,inf)>e    
    D=max(diag(A*x0-b),0);
    D(D>0)=1;
    d=x0./(A'*D*A*x0+det);
    p1=-d.*g;
    c=(g'*g)/(g'*A'*A*g);
    if x0-c*g>0
        c=(g'*g)/(g'*A'*A*g);
    else
        C=x0./g;
        c=min(C);
    end
    p2=-c*g;
    p=p1+p2;
   
   t=ffindt(A,p1,p2,g,D);
   s=t*(p2-p1)+p1;
    x0=x0+s;
    g=fdetq(A,b,x0);  
    i=i+1;
    if i==100        
        break;
    end
end
x=x0;
end

