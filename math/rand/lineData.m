% generate point in the range x=[xmin xmax] y=[ymin ymax]
% b-a(1)x = A(2)y A must two column
function [d,out]=lineData(A,b,x,y)
[m,n]=size(A);
d=repmat([x(1),y(1),x(2),y(2)],[m,1]);
for i=1:m
    if abs(A(i,1))>1e-8 && abs(A(i,2))>1e-8
        % x1 = (b - a(2)ymin)/a(1) 
        d(i,1)=(b(i,1)-A(i,2)*y(1))/A(i,1);
        % x1 = (b - a(2)ymax)/a(1) 
        d(i,2)=(b(i,1)-A(i,2)*y(2))/A(i,1);
        % y1 = (b - a(1)xmix)/a(2) 
        d(i,3)=(b(i,1)-A(i,1)*x(1))/A(i,2);
        % y1 = (b - a(1)xmix)/a(2) 
        d(i,4)=(b(i,1)-A(i,1)*x(2))/A(i,2);
        %     , ymin;
        %     ,ymax;
        %xmin,     ;
        %xmax,
        dd=[d(i,1),y(1);d(i,2),y(2);x(1),d(i,3);x(2),d(i,4)];
        % take the middle point to satisfy range 
        dd=sortrows(dd,1);
        % x 
        d(i,1)=dd(2,1);
        d(i,2)=dd(3,1);
        % y 
        d(i,3)=dd(2,2);
        d(i,4)=dd(3,2);
        % d= [x1,x2,y1,y2;] every rows for every line
        % line horizon
    elseif abs(A(i,1))>1e-8
        d(i,1)=b(i,1)/A(i,1);
        d(i,2)=b(i,1)/A(i,1);
        d(i,3)=y(1);
        d(i,4)=y(2);   
    elseif abs(A(i,2))>1e-8
        d(i,1)=x(1);
        d(i,2)=x(2);
        d(i,3)=b(i,1)/A(i,2);
        d(i,4)=b(i,1)/A(i,2);   
    end  
    out = 
end

% valid data
% A=[1,1;-1,-1;1,0;6,3];
% b=[1;1;-0.5;-2];
% % output
% p=[1,0,0,1;-1,0,0,-1;-1/2,0,-1/2,1;-1/2,0,0,-2/3];
% d=lineData(A,b);
% d-p



