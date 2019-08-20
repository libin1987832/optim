function [real,se,both]=diffActive(wset,realA,sep)
both=[];
real=[];
se=[];
for q=1:size(wset,2)
    if sum(find(ismember(realA,wset(q))>0))>0 && sum(find(ismember(sep,wset(q))>0))>0
        fprintf('%dV',wset(q));
        both=[both wset(q)];
    else if sum(find(ismember(realA,wset(q))>0))>0
            fprintf('%dA',wset(q));
            real=[real wset(q)];
        else if sum(find(ismember(sep,wset(q))>0))>0
                fprintf('%dS',wset(q));
                se=[se wset(q)];
            else
                fprintf('%d',wset(q));
            end
        end
        
    end
    if q~=size(wset,2)
        fprintf(' ');
    end
end