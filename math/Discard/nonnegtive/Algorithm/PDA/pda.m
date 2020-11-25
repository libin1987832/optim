function x0=pda(A,b,n)
x0=100*ones(n,1);
active=[];
while(1)
    [x1,r1,active,num]=nextIter(A,b,active,n,x0);
    x0=x1;
    disp(['active:',num2str(active),'num:',num2str(num)]);
    z=x1.*r1;
    if norm(z)<0.0001
        break;
    end
end


