%����������������
function [x,fk]=fsearchx(A,b,x0,e,det,maxit)
%ʹ�ò�����������x
g=fdetq(A,b,x0);
index = 1;
%while norm(x0.*g,inf)>e
while norm(x0.*g,inf)>e||min(g)<-e
    index =index +1;
    if index >maxit
        break;
    end
    D=max(diag(A*x0-b),0);
    D(D>0)=1;
		% d=x/(A'DAx+det)
    d=x0./(A'*D*A*x0+det);
		%p=-d.*det
    p=-d.*g;

    %a=fsearcha(A,b,x0,p); %ʹ��ԭ������������
    %a=fwolfepowersearcha(A,b,x0,p);%ʹ��wolfepower��������������������
%    [alpha,x0,fx0,g]=wolfe(A, b, x0, p, 0.9); 
    [alpha, minf, knot] = arraySpiecewise(A,b,x0,p);  
     x0 = x0 + alpha*p; 
end
x=x0;
fk=fq(A, b,x);
end

