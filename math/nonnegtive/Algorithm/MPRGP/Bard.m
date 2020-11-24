function  [w,z,retcode] = Bard(M,q,pivtol,maxits)

if nargin<3, pivtol = 1e-8; maxits = 1e4; end;
if nargin<4, maxits = 1e3; end;
n = length(q);
if size(M)~=[n n]; error('Matrices are not compatible'); end;

rayTerm = false;
loopcount = 0;
if min(q)>=0 % If all elements are positive a trivial solution exists
    %  As w - Mz = q, if q >= 0 then w = q a, and z = 0
    w = q;
    z = zeros(size(q));
else
    dimen = size(M,1); % Number of rows
    I = eye(dimen);
    % Let artificial variable enter the basis
    basis = 1:dimen; % A set of row indices in the tableau
    Mq = [M,q];
    % Perform complementary pivoting
    while min(q)<0 && loopcount < maxits
        [~,locat] = min(q); % Row of minimum element in column
        % 2*dimen+1 (last of tableau)

        P = -Mq(:,locat)/Mq(locat,locat);
        P(locat) = -1/Mq(locat,locat);
        Mq(:,locat) = -I(:,locat);
        
        pivot = Mq(locat,:);
        Mq = Mq + P*pivot;
        Mq(locat,:) = P(locat)*pivot;
        oldVar = basis(locat);
        % Select next candidate for entering the basis
        if oldVar > dimen
            basis(locat) = oldVar - dimen;
        else
            basis(locat) = oldVar + dimen;
        end
    end
    % Return the solution to the LCP
    vars = zeros(2*dimen+1,1);
    vars(basis) = tableau(:,dimen+2).';
    w = vars(1:dimen,1);
    z = vars(dimen+1:2*dimen,1);
end

if rayTerm
    retcode = [2, loopcount];  % Ray termination
else
    retcode = [1, loopcount];  % Success
end
