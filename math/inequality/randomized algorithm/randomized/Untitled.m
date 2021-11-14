        alphSort=sort([minrcol;plaha(minrcol < plaha);malpha(malpha < maxrcol);maxrcol]);
       for j = 1:size(alphSort,1)
           rj=b-alphSort(j)*row;
           df=row'*rj;
           if df >0
                break;
           end
       end
       
              Rrepmat = repmat(b,1,size(alphSort,1)) - row*alphSort';
       Rrepmat(Rrepmat<0)=0
       
       
       %        accountmin = minrcol < plaha;
%        active1 = [plaha(accountmin),-ones(sum(accountmin),1)];
%        accountmax = malpha < maxrcol;
%        active2 = [malpha(accountmax),ones(sum(accountmax),1)];
  
        alphSort=sort([minrcol;plaha(minrcol < plaha);malpha(malpha < maxrcol);maxrcol]);
       for j = 1:size(alphSort,1)
           rj=b-alphSort(j)*row;
           df=row'*rj;
           if df >0
                break;
           end
       end
     %  active = sortrows([active1,active2],'ascend');
       AN_r+r(i)*A(i,:)' - active*A(:,i)'*col