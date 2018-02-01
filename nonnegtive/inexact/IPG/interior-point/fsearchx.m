%靠靠靠靠靠靠靠靠
function [x,fk]=fsearchx(A,b,x0,e,det)
%使用步长来找最优x
g=fdetq(A,b,x0);
%while norm(x0.*g,inf)>e
while norm(x0.*g,inf)>e||min(g)<-e

    D=max(diag(A*x0-b),0);
    D(D>0)=1;
		% d=x/(A'DAx+det)
    d=x0./(A'*D*A*x0+det);
		%p=-d.*det
    p=-d.*g;
    %a=fsearcha(A,b,x0,p); %使用原步长搜索方法
    %a=fwolfepowersearcha(A,b,x0,p);%使用wolfepower线型搜索方法搜索步长
    [alpha,x0,fx0,g]=wolfe(A, b, x0, p, 0.9); 
end
x=x0;
fk=fq(A, b,x);
end

