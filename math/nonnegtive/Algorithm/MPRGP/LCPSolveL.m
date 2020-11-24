function [w,z,retcode] = LCPSolveL(M,q,pivtol,maxits)

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
    
    % Let artificial variable enter the basis
    basis = 1:dimen; % A set of row indices in the tableau
    nonbasis = (dimen+1):(2*dimen+1);
    % 选择最小元素，作为输入变量
    [~,locat] = min(q); % Row of minimum element in column 
                                     % 2*dimen+1 (last of tableau)
    Ie=zeros(dimen,1);
    Ie(locat) = -1;
   % Create initial tableau
    % |   |    |    |   |
    % | I | M | e | q |
    tableau = [M, Ie, q];
    P = -ones(dimen,1);
    % 选择基向量的索引                               
    basis(locat) = 2*dimen+1; % Replace that index with the column
    nonbasis(dimen+1) = locat;
    % 选择基向量的索引   
    cand = locat + dimen;
    
    pivot =  tableau(locat,:);
     
    % From each column subtract the column 2*dimen+1, multiplied by pivot
    % P = I - pivot
    tableau = tableau + P*pivot; 
    tableau(locat,:) = P(locat)*pivot;
    % Perform complementary pivoting
    while max(basis) == 2*dimen+1 && loopcount < maxits
        cand3 = find(nonbasis == cand);
        loopcount = loopcount + 1;
        eMs = tableau(:,cand3); % This is used to check convergence (pivtol)
        missmask = eMs >= 0;  % Check if elements of eMs are less than zero
        quots = tableau(:,dimen+2)./eMs;
        quots(missmask) = -Inf;
        [~,locat] = max(quots);
        % Check if at least one element is not missing
        if  sum(missmask)~=dimen && abs(eMs(locat)) > pivtol 
            P = -tableau(:,cand3)/tableau(locat,cand3);
            P(locat) = -1/tableau(locat,cand3);
            tableau(:,cand3) = zeros(dimen,1); 
            tableau(locat,cand3) = -1; 
            % Reduce tableau
            pivot =  tableau(locat,:);
            tableau = tableau + P*pivot;
            tableau(locat,:) = P(locat)*pivot;
            
            oldVar = basis(locat);
            % New variable enters the basis
            basis(locat) = cand;
            nonbasis(cand3) = oldVar;
            % Select next candidate for entering the basis
            if oldVar > dimen
                cand = oldVar - dimen;
            else
                cand = oldVar + dimen;
            end
        else
            rayTerm = true; % Problem was solved
            break % Break out of the while loop
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

