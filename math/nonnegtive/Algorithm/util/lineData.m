function d=lineData(A,b,x,y)
[m,n]=size(A);
d=repmat([x(1),y(1),x(2),y(2)],[m,1]);
for i=1:m
    if abs(A(i,1))>1e-8
        d(i,1)=(b(i,1)-A(i,2)*y(1))/A(i,1);
    else
        d(i,1)=1;
        if abs(A(i,2))>1e-8
            d(i,2)=(b(i,1)-A(i,1)*x(1))/A(i,2);
        else
            print('zores')
        end
    end
    if abs(A(i,2))>1e-8
        d(i,4)=(b(i,1)-A(i,1)*x(2))/A(i,2);
    else
        d(i,4)=1;
        if abs(A(i,1))>1e-8
            d(i,3)=(b(i,1)-A(i,2)*y(1))/A(i,1);
        else
            print('zores')
        end
    end
end

% valid data
% A=[1,1;-1,-1;1,0;6,3];
% b=[1;1;-0.5;-2];
% % output
% p=[1,0,0,1;-1,0,0,-1;-1/2,0,-1/2,1;-1/2,0,0,-2/3];
% d=lineData(A,b);
% d-p



