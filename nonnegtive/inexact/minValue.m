% 真真真真真真真真真真真真真真真真�
% x 真真真  scale 真真真真真真 number 真真真 ||b-Ax|| 真真
function [minNum,v]=minValue(x,scale,number,A,b,f,c)
[m,n]=size(x);
fx=f(A,b,x);
v=zeros(m,number);
minNum=0;
index=0;
while 1
    %昧字恢伏方忖
		% 真真真真�
    xi=x+(2*rand(m,n)-ones(m,n))*scale;
		% 真真真真真真
    g=c(xi);
    %諾怎埃崩訳周
		% 真
    if g>0
				% 真真真		
        fxi=f(A,b,xi);
        %曳公協恷弌峙珊勣�a�
				% 真真真真真真� 真真真
        if fxi<fx
            minNum=minNum+1;
            v(:,minNum)=xi;
        end
        %峪勣憲栽埃崩議峙嘉嬬麻嗤丼峙
        index=index+1;
        %糞刮峙曳公協議勣謹 祥唯峭登僅
        if index>number
            break;
        end
    end
end; 
