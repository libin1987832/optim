% interv [0,1] for choose a the projectedsearch but less one 
function [alpha,x0,knot, retcode] = arraySpiece(A,b,x0,p,tol,maxits)
display = true;
if nargin<5, tol = 1e-15; maxits = 1e4; end
[m,n] = size(A);
nonzerou = p < 0;
knot = -x0(nonzerou)./p(nonzerou);
knot = sort(knot( knot > tol & knot < 1));
knot = [0;knot;1];
knot = unique(knot);
left = 1;
right = length(knot);
loopcount = 1;
alpha = 1;
retcode = [1, 0];
%for i = 2:length(aranges)
while left+1  <= right && loopcount < maxits
    loopcount = loopcount + 1;
    % floor down so that i+1
    i = floor(0.5*(left + right));
    % take the middle value for active set
    x = x0 + 0.5*(knot(i) + knot(i+1))*p;
    % x active setss
    Iu = x < 1e-15;
    [alpha,retcode] = spiecewise(A(:,~Iu),b,p(~Iu),x0(~Iu),knot( i ) , knot( i+1 ));
    if retcode(1) == 1 && retcode(2) == 0
%         if display
%             fprintf("newton arrayspiece min alpha in the left, so stop\n");
%         end
        break;
    end
    if retcode(1) == 2 && retcode(2) == 0
%         if display
%             fprintf("newton arrayspiece p is close to zeros\n");
%         end
        break;
    end
    if retcode(1) == 1 && retcode(2) ==1
        right = i;
    end
    if retcode(1) == 2 && retcode(2) ==1
        left = i;
    end
end
xk = x0 + alpha * p;
xk(xk<0) =0;
if display
    knotL = length(knot);
    [rpk0, normr0, xmin0, Ar0, normKKT0 , face11, face12] = kktResidual(A, b, x0, [], 1);
    [rpk1, normr1, xmin1, Ar1, normKKT1 , face21, face22] = kktResidual(A, b, xk , [], 1);
    fprintf('___projected(0,1): alpha(%g),knot_numb(%d),normB(%g),normF(%g),gB(%g),gF(%g),xa(%d,%d),ba(%d,%d)',...
        alpha,knotL,normr0,normr1,Ar0,Ar1,face11, face21,face12, face22);
end
%minf = funmin(A,b,x0,p,alpha);