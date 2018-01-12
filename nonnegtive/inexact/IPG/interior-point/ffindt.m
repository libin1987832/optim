function t=ffindt(A,p1,p2,g,D)
%b=0.9;
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%(1/2)*t^2*(p2-p1)'*A'*A*(p2-p1)+t*p1'*A'*A*(p2-p1)+t*g'*(p2-p1)+(1/2)*p1'*A'*D*A*p1+p1'*g-b*((1/2)*p2'*A'*D*A*p2+p2'*g)
t1=(1/2)*(p2-p1)'*A'*A*(p2-p1);%a

t2=p1'*A'*A*(p2-p1)+g'*(p2-p1);%b

t3=(1/2)*p1'*A'*D*A*p1+p1'*g-0.9*((1/2)*p2'*A'*D*A*p2+p2'*g);%c

%t=-1/2*(t2-(t2^2-4*t1*t3)^(1/2))/t1;
t=-1/2*(t2+(t2^2-4*t1*t3)^(1/2))/t1;
end

