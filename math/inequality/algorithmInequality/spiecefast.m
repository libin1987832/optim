function steplength=spiecefast(A,b,p,x0)
if p'*p<1e-15
    steplength = 0;
    return
end
r=b-A*x0;
ApIu = A*p;
knotr=r./ApIu;
knotr = sort( knotr( knotr >=0 & knotr <= 1));
knotr= unique([0; knotr ; 1]);
leftb = 1;
rightb = length(knotr);
loopcountB = 0;
maxits = 300;
while leftb + 1  <= rightb && loopcountB < maxits
    loopcountB = loopcountB + 1;
    if loopcountB == 1
        knotri = rightb - 1;
    elseif loopcountB == 2
        knotri = leftb;
    else
        knotri = floor(0.5*(leftb + rightb));
    end
    
    alpha = knotr(knotri);
    beta = knotr(knotri + 1);
    
    %                 testalphab(loopcountB) = alpha;
    
    % middle point give the active set
    rknot = r  - 0.5 * (knotr(knotri)+ knotr(knotri + 1)) * ApIu;
    % the derive Ap'(r-alpha*Ap)=0
    Apr = ApIu(rknot > 0);
    Ar = Apr' * r(rknot > 0);
    % alpha = (Ap*r)^T/(Ap'*Ap)
    steplength = Ar / ( Apr' * Apr );
    % the right point
    if steplength >= beta
        steplength = beta;
        if loopcountB == 1
%             retcode = [1,1];
            break;
        else
            leftb = knotri+1;
%             retcode = [1,0];
        end
    elseif steplength <= alpha
        steplength = alpha;
        rightb = knotri;
        if loopcountB == 2
%             retcode = [1,2];
            break;
        end
    else
        break;
    end
end