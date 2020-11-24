function [w,z,retcode] = Graves(M,q,pivtol,maxits)
if nargin<3, pivtol = 1e-8; maxits = 1e4; end;
if nargin<4, maxits = 1e3; end;
n = length(q);
if size(M)~=[n n]; error('Matrices are not compatible'); end;

loopcount = 1;

if min(q)>=0 % If all elements are positive a trivial solution exists
    %  As w - Mz = q, if q >= 0 then w = q a, and z = 0
    w = q;
    z = zeros(size(q));
     retcode = [1, loopcount]; 
else
    dimen = size(M,1); % Number of rows
    I = eye(dimen);
    % Let artificial variable enter the basis
    basis = 1:dimen; % A set of row indices in the tableau
    MqI = [M,I,q];
    % Perform complementary pivoting
    while min(MqI(:,end))<0 && loopcount < maxits
        loopcount = loopcount + 1;
        indexq = MqI(:,end) < 0;
      
        if sum(indexq) ==1
            locat = find(indexq>0);
        else
            tMqI = MqI(:,dimen+1:2*dimen);
            tMqI(~indexq,:) = -Inf;
            test = tMqI(indexq,1:dimen)./repmat(MqI(indexq,end),1,dimen);
            y = any(diff(test'),2);
            index = find(y ~= 0);
            [~,locat] = max(test(:,index(1)));
        end
        % 2*dimen+1 (last of tableau)

        P = -MqI(:,locat)/MqI(locat,locat);
        P(locat) = -1/MqI(locat,locat);
        MqI(:,locat) = -I(:,locat);
        
        pivot = MqI(locat,:);
        MqI = MqI + P*pivot;
        MqI(locat,:) = P(locat)*pivot;
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
    vars(basis) = MqI(:,end).';
    w = vars(1:dimen,1);
    z = vars(dimen+1:2*dimen,1);
    retcode = [1, loopcount]; 
end

% 
% if rayTerm
%     retcode = [2, loopcount];  % Ray termination
% else
%     retcode = [1, loopcount];  % Success
% end
