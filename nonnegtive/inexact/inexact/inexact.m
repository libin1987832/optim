% ¿¿¿¿¿
% ||(b-Ax)+|| x0¿¿ belt¿0.5¿ ¿¿¿¿¿¿¿¿¿¿¿¿¿¿ u¿0-1¿ ¿¿¿¿¿¿¿¿¿
% x(k+1)=(xk+delt^m*a*s) ||bk-Ax(k+1)||<=||y||-2u*s*(x(k+1)-xk)
function [x0,f1]=inexact(b,A,x0,belt,u)
[msize,nsize]=size(x0);
index=0;
while 1
		% ¿¿¿¿¿¿¿¿		
    y0=b-A*x0;
    y0(y0<=0)=0;
    % ¿¿¿¿¿¿¿¿¿¿
    b0=y0+A*x0;
    s0=A'*y0;
    % ¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿
    sN=s0;
    alph=sN'*sN/((A*sN)'*(A*sN));
    mk=0;
    lamd=diag(ones(msize,1));
    for i=1:msize
        if 0<x0(i) & s0(i)<0
            lamd(i,i)=-1*x0(i)/s0(i);
        end
    end
%     %²âÊÔÊÇ·ñ´æÔÚÕâÑùµÄÏµÊý
%     xxxx=[];
%     yyyy=[];
%     for tt=1:1000
%         x1=x0+1/100000*tt*s0;
% %         x1(x1<0)=0;
%         left=(b0-A*x1)'*(b0-A*x1);
%         right=y0'*y0-2*u*s0'*(x1-x0);
%         ff=left-right;
%         xxxx=[xxxx,1/100000*tt];
%         yyyy=[yyyy,ff];
%     end
%       plot(xxxx,yyyy);
%       %%%%%%
    % ¿¿¿¿¿¿¿¿¿¿¿¿x(k+1)=(xk+delt^m*a*s) ||bk-Ax(k+1)||<=||y||-2u*s*(x(k+1)-xk)¿¿
    while 1
				% ¿¿¿belt,belt^2,belt^3¿¿¿¿¿¿		
        x1=x0+belt^mk*alph*s0;
%         x1=x0+belt^mk*lamd*s0;
        x1(x1<0)=0;
        left=(b0-A*x1)'*(b0-A*x1);
        right=y0'*y0-2*u*s0'*(x1-x0);
				% ¿¿¿¿¿
        if left<=right
            break;
        end
        mk=mk+1;
    end
		% ¿¿¿¿¿¿¿
    f0=b-A*x0;
    f0(f0<0)=0;
    f0=0.5*(f0'*f0);
		% ¿¿¿¿¿¿¿¿
    f1=b-A*x1;
    f1(f1<0)=0;
    f1=0.5*(f1'*f1);
    fprintf('index:%d,The mk is %f; x0(%f,%f),x1(%f,%f); f0,f1:%f,%f\n',index,mk,x0,x1,f0,f1);
    index=index+1;
    x0=x1;
    % ¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿
    if norm(f1-f0)<0.00001
        break;
    end
end
