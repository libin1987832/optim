% ¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿
% x ¿¿¿¿¿¿  scale ¿¿¿¿¿¿¿¿¿¿¿¿ number ¿¿¿¿¿¿ ||b-Ax|| ¿¿¿¿
function [minNum,v]=minValue(x,scale,number,A,b,f,c)
[m,n]=size(x);
fx=f(A,b,x);
v=zeros(m,number);
minNum=0;
index=0;
while 1
    %Ëæ»ú²úÉúÊý×Ö
		% ¿¿¿¿¿¿¿¿¿
    xi=x+(2*rand(m,n)-ones(m,n))*scale;
		% ¿¿¿¿¿¿¿¿¿¿¿¿
    g=c(xi);
    %Âú×ãÔ¼ÊøÌõ¼þ
		% ¿¿
    if g>0
				% ¿¿¿¿¿¿		
        fxi=f(A,b,xi);
        %±È¸ø¶¨×îÐ¡Öµ»¹ÒªÐa¡
				% ¿¿¿¿¿¿¿¿¿¿¿¿¿ ¿¿¿¿¿¿
        if fxi<fx
            minNum=minNum+1;
            v(:,minNum)=xi;
        end
        %Ö»Òª·ûºÏÔ¼ÊøµÄÖµ²ÅÄÜËãÓÐÐ§Öµ
        index=index+1;
        %ÊµÑéÖµ±È¸ø¶¨µÄÒª¶à ¾ÍÍ£Ö¹ÅÐ¶Ï
        if index>number
            break;
        end
    end
end; 
