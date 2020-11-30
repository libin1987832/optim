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
    f1 = funmin(A,b,x0,p,aranges( i - 1 ));
    f2 = funmin(A,b,x0,p,aranges( i ));
    if f1 < f2
        alpha = aranges(i-1);
        minf = f1;
        return;
    end
    % take the middle value for active set
    u = p;
    x = x0 + aranges( i - 1 )*u;
    % x active setss
    Iu = x < 1e-15;
    x(Iu) = 0;
    u(Iu) = 0;
    [alpha,retcode] = spiecewise(A,b,u,x,aranges( i )-aranges( i - 1 ));
    if retcode(1) == 1 && retcode(2) == 0
        break;
    end
end
minf = funmin(A,b,x0,p,alpha);