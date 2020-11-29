%���������������� [xk,resvec,arvec,face1vec,face2vec,tf]
function [xk,resvec,arvec,face1vec,face2vec,tf] = fsearchx(A,b,x0,e,det,maxit)
%ʹ�ò�����������x
t=clock;
g=fdetq(A,b,x0);
index = 1;

 [rpk, normr, xmin, Ar, normKKT, face1, face2] = kktResidual(A, b, x0 , [], 1);
    
% the residual vector
resvec = zeros(1,maxit);
% the normal gradient 
arvec = zeros(1,maxit);
% subspace minization
itersm = zeros(1,maxit);
% face b-Ax
face1vec = zeros(1,maxit);
% face x
face2vec = zeros(1,maxit);
resvec(1) = xmin;
arvec(1) = normKKT;
face1vec(1) =face1;
face2vec(1) =face2;

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
     
     [rpk, normr, xmin, Ar, normKKT, face1, face2] = kktResidual(A, b, x0 , [], 1);
    % record the value of objection function
    resvec(index) = normr;
    % record the value of the gradient function
    arvec(index) = normKKT;
    face1vec(index) = face1;
    face2vec(index) = face2;
end
xk=x0;
tf = etime(clock,t);
end

