% interv [0,1] for choose a 
function [alpha, knot, retcode] = arraySpiece(A,b,x0,p,tol,maxits)
display = true;
if nargin<5, tol = 1e-15; maxits = 1e4; end
[m,n] = size(A);
nonzerou = p < 0;
knot = -x0(nonzerou)./p(nonzerou);
knot = sort(knot( knot > tol & knot < 1));
knot = [0;knot;1];
knot = unique(knot);
left = 1;
right = length(knot)-1;
loopcount = 1;
alpha = 1;
retcode = [1, 0];
%for i = 2:length(aranges)
while left+1  < right && loopcount < maxits
    loopcount = loopcount + 1;
    i = floor(0.5*(left + right));
    % take the middle value for active set
    x = x0 + 0.5*(knot(i - 1) + knot(i))*p;
    % x active setss
    Iu = x < 1e-15;
    [alpha,retcode] = spiecewise(A(:,~Iu),b,p(~Iu),x0(~Iu),knot( i-1 ) , knot( i ));
    if retcode(1) == 1 && retcode(2) == 0
        if display
            fprintf("min alpha in the left, so stop\n");
        end
        break;
    end
    if retcode(1) == 2 && retcode(2) == 0
        if display
            fprintf("p is close to zeros\n");
        end
        break;
    end
    if retcode(1) == 1 && retcode(2) ==1
        right = i;
    end
    if retcode(1) == 2 && retcode(2) ==1
        left = i;
    end
end
%minf = funmin(A,b,x0,p,alpha);