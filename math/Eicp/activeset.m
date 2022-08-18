function [act,nact] = activeset(x,y,c)
[m,n]=size(x);
% for i=1:m
%    if(c*x(i)+y(i)>0)
%        act(i)=1;
%        nact(i)=0;
%    else
%        act(i)=0;
%        nact(i)=1;
%    end 
% end

for i=1:m
   if(x(i)==0&&y(i)>0)
       act(i)=1;
       nact(i)=0;
   else
       act(i)=0;
       nact(i)=1;
   end 
end