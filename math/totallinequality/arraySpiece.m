% interv [0,1] for choose a 
function [alpha, minf, aranges, retcode] = arraySpiece(A,b,x0,p,tol,maxits)
if nargin<5, tol = 1e-15; maxits = 1e4; end
[m,n] = size(A);
nonzerou = p < 0;
arange = -x0(nonzerou)./p(nonzerou);
aranges = sort(arange(arange > tol & arange < 1));
aranges = [0;aranges;1];
aranges = unique(aranges);
for i = 2:length(aranges)
    % take the middle value for active set
    x = x0 + 0.5*(aranges(i - 1)+aranges(i))*p;
    % x active setss
    Iu = x < 1e-15;
    [alpha,retcode] = spiecewise(A(:,~Iu),b,p(~Iu),x0(~Iu),aranges( i-1 ) , aranges( i ));
    if retcode(1) == 1 
        break;
    end
end
minf = funmin(A,b,x0,p,alpha);