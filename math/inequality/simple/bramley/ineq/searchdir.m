%
%        function [d,R,error,E] = searchdir(A,I,Iold,R,r,method,E)
%
%  searchdir uses one of the 4 methods (SVD if method == 1;
%  QR with pivoting if method == 2; or up/downdate if method ==3)
%  to find the search direction
%
% E is the permutation vector (if any) from a previous QR 
%  factorization.
% error = 1 means that update or downdate did not work
%
%   Beata Winnicka
%   Dept of Computer Science
%   Indiana University - Bloomington
%   1 October 1993 
%
% Modified: R. Bramley
%           Thu Mar 31 16:03:17 EST 1994
%           To use space saver version of QR factorization.
%-------------------------------------------------------------------------------

      function [d,R,error,E] = searchdir(A,I,Iold,r,method,E)
      [m,n] = size(A(I,:));
	  error = 0;
%
%  Method = 0 means to use complete orthogonal factorization
%
      if (method == 0) 

        [Q,R,E] = qr(A(I,:),0);
%  Find numerical rank of R:
        for p = 1:min(n,m);
            if abs(R(p,p)) < abs(R(1,1))*max(size(A(I,:))*eps);
               disp(['                SD info: Rank deficiency detected; p = ',num2str(p-1)])
               p = p - 1;
               break;
            end;
        end;
        d = Q'*r(I);
        
        H = zeros(n-p+1,p);

        for i = p:-1:1
           h = house(R(i,[i,p+1:n]));
           R(1:i,[i,p+1:n]) = R(1:i,[i,p+1:n]) - (R(1:i,[i,p+1:n])*h)*h';
           H(:,i) = h;
        end  %  for i = p:-1:1

        d = R(1:p,1:p)\d(1:p);
        d = [d;zeros(n-p,1)];

        for i = 1:p
            d([i,p+1:n]) = d([i,p+1:n]) - (H(:,i)'*d([i,p+1:n]))*H(:,i);
        end  %  for i = 1:p 

        d(E) = d;
    
      end  %  if (method == 0) 
    
%
%   Method = 1 means use singular value decomposition
%
      if (method == 1) 
              [u,s,v] = svd(A(I,:));
              s = diag(s);
              utr = u'*r(I);
              y = zeros(n,1);
              [n1,n2] = size(I);
              for kk = 1:min(n,n1);
                    if s(kk) > s(1)*max(size(A(I,:))*eps);
                      y(kk) = utr(kk)/s(kk);
                    else
                      disp(['                SD info: Rank deficiency detected; rank = ',num2str(kk-1)])
                      break;
                    end;
              end;
              d = v*y;
              R=0;
      end 
%
%  Method = 2 means to use QR with column pivoting 
%
    if (method == 2)
        [Q,R,E] = qr(A(I,:),0);
%  Find numerical rank of R:
        for p = 1:min(n,m);
            if abs(R(p,p)) < abs(R(1,1))*max(size(A(I,:))*eps);
               disp(['                SD info: Rank deficiency detected; p = ',num2str(p-1)])
               p = p - 1;
               break;
            end;
        end;
        d = Q'*r(I);
        d = R(1:p,1:p)\d(1:p);
        d = [d;zeros(n-p,1)];
        d(E) = d;
    end 
%
%  Method = 3 means use Cholesky factorization update/downdate
    if (method == 3)
       Ioldnew=find(Iold >0);
%  Check to make sure there are as many active rows as columns; otherwise
%  R is singular and update/downdate may not work.
       if(length(I) < n),
          error = 1; return;
       end
       len=max([I;Ioldnew]);
       h1=zeros(1,len);
       h2=zeros(1,len);
       h1(I)=ones(size(I));
       h2(Ioldnew)=ones(size(Ioldnew));
       Idown=find((h1-h2)<0);
       Iup=find((h2-h1)<0);
       disp(['                SD info: need to downdate ',num2str(length(Idown)),' rows'])
       disp(['                SD info: need to update ',num2str(length(Iup)),' rows'])

%    check if num flops for up/downdate better than regular QR
%
       if((4*n*n+10*n)*(length(Iup)+length(Idown)) >= m*n*n - n*n*n/3)
          error = 1;
          disp(['                SD info: up/downdate flops larger than regular QR'])
          return
       end
%
       for jup=1:length(Iup)
%           [error, R] = update(R, A(Iup(jup),:)*E);
           temp = A(Iup(jup),E);
           [error, R] = update(R, temp);
           if error == 1
             disp(['                SD info: update failed '])
             return
           end % error == 1
       end % for jup=1:length(Iup)
%
       for jdown=1:length(Idown)
%           [error, R] = downdate(R, A(Idown(jdown),:)*E);
           temp = A(Idown(jdown),E);
           [error, R] = downdate(R, temp);
           if error == 1  
              disp(['                SD info: downdate failed '])
              return
           end % error == 1
       end % for jdown=1:length(Idown)
%
%   calculate new d
%
        [p1,p2]=size(R);
        p = min(p1,p2);
%        d = E*(R(1:p,1:p)\(R(1:p,1:p)'\(E'*(A(I,:)'*r(I)))));
        temp = R(1:p,1:p)\(R(1:p,1:p)'\(A(I,E)'*r(I)));
        d(E) = temp;
        d = d(:);
   end % if (method == 3)

   return
