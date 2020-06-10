function xkA= splitS_ourres(Q,d,s,x0,iter)
x1=x0;
x2=x0;
x3=x0;
if iter<1
    xk=x0;
    xkA=[xkA xk];
else
    [m,n]=size(Q);
    tol=1e-12;
    % x0=zeros(n,1);
    for i=1:iter
        %xk=zeros(n,1);
        x=x0;
        for j=1:n
            old_xi = x(j);
            ri     = d(j) + Q(j,:)*x;
            Aii    = Q(j,j);           
            x(j) = max( 0, old_xi - (ri / Aii) );
        end

%         for j=1:n
%             sd=s*(d(j)+Q(j,:)*[xk(1:(j-1));x0(j:n)])/Q(j,j);
%             if x0(j)>sd
%                 xk(j)=x0(j)-sd;
%             end
%         end
%         xn1=x0;
%         x0=xk;
        x0=x;
        x1=x2;
        x2=x3;
        x3=x0;
    end
end
xkA=[x1 x2 x3];