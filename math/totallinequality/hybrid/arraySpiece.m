% interv [0,maxknot] for choose a the projectedsearch but less maxknot
function [alpha,xk,knot, retcode] = arraySpiece(A,b,x0,p,tol,maxits)
display = false;
if nargin<5, tol = 1e-15; maxits = 1e4; end
[m,n] = size(A);
nonzerou = p < 0;
knot = -x0(nonzerou)./p(nonzerou);
maxknot = 1e5;
knot = sort(knot( knot > tol & knot < 1e5));
knot = [0;knot;1e5];
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
%     if alpha <1e-10;
%         xa = [1e-12:1e-14:2*1e-12];
%         ya = arrayfun(@(alpha) funmin(A,b,x0,p,alpha), xa);
%         pxy={};
%         pxy(1).X = xa;
%         pxy(1).Y = ya;
%         figure
%         hold on
%         p1 = arrayfun(@(a) plot(a.X,a.Y),pxy);
%     end
    knotL = length(knot);
    [rpk0, normr0, xmin0, Ar0, normKKT0 , face11, face12] = kktResidual(A, b, x0, [], 1);
    [rpk1, normr1, xmin1, Ar1, normKKT1 , face21, face22] = kktResidual(A, b, xk , [], 1);
    fprintf('___projected(0,1): alpha_num(%g,%d),normBF(%g,%g),gBFP(%g,%g,%g),xa(%d,%d,%d,%d)\n',...
        alpha,knotL,normr0,normr1,norm(Ar0),norm(Ar1),Ar0'*p/(norm(Ar0)*norm(p)),face11, face21,face12, face22);
end
%minf = funmin(A,b,x0,p,alpha);