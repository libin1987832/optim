%  function[error, R] = downdate(R,x)
%
%   Performs single row downdate of R
%
%   Beata Winnicka
%   Dept of Computer Science
%   Indiana University - Bloomington
%   1 October 1993 
%
%
   function[error, R] = downdate(R,x)

   eps=1.0e-5;
   [dum,p]=size(R);
   s=zeros(p,1);
   c=zeros(1,p);
%   flopsnow = flops;
   s=R(1:p, 1:p)'\x';
%   flopsnow=flops - flopsnow;
   snorm = norm(s);
   if(snorm>1-eps)
        error=1;
        return;
   end
   alpha=sqrt(1.0-snorm*snorm);

   for ii=1:p
       i=p-ii+1;
       scale=alpha+abs(s(i));
       a=alpha/scale;
       b=s(i)/scale;
       tnorm=sqrt(a*a+b*b);
       c(i)=a/tnorm;
       s(i)=b/tnorm;
       alpha=scale*tnorm;
   end
   for j=1:p
       xx=0.0;
       for ii=1:j
           i=j-ii+1;
           t=c(i)*xx+s(i)*R(i,j);
           R(i,j)=c(i)*R(i,j)-s(i)*xx;
           xx=t;
        end
   end
   error=0;
