% 判断是否所有的分量为0
function g=constrain(x)
     a=find(x<0);
     g=1;
     if sum(a)>0
         g=0;
     end
end
     
