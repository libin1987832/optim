% interval [0,1] for choosing steps but discarded, since the algorithm is the sequence
function [alpha, minf, knot,retcode] = arraySpiecewise(A,b,x0,p,tol,maxits)
if nargin<5, tol = 1e-15; maxits = 1e4; end
[m,n] = size(A);
% find the valid component
nonzerou = p < 0;
arange = -x0(nonzerou)./p(nonzerou);
% sort
aranges = sort(arange(arange > tol & arange < 1));
aranges = [0;aranges;1];
aranges = unique(aranges);
all = [];
% the sequence for searching
for i = 2:length(aranges)
    % take the middle value for active set
    u = p;
    x = x0 + aranges( i - 1 )*u;
    % x active setss
    Iu = x < 1e-15;
    x(Iu) = 0;
    % z active set
    r = b - A * x;
    u(Iu) = 0;
    % (b-Ax-Au)
    ap = A*u;
    % find the valid component for the residual
    ai=r./ap;
    alength = aranges( i ) - aranges( i -1);
    % we are concentrated on the valid interval 
    as=sort(ai(ai > 0 & alength >ai));
    % displace the distance
    ass = aranges(i-1) + as;
    % all knots
    ass = [aranges( i-1 ) ass' aranges(i)]; 
    % just for debug
    all = [all,ass];
end
% summary all possible value
knot = unique(all);
% check the interval in turn
for i = 2 : length(knot)
    % the object function for i-1
    f1 = funmin(A,b,x0,p,knot( i - 1 ));
    % the object function for i
    f2 = funmin(A,b,x0,p,knot( i ));
    % It means the objection function is increasing, so stop
    if f1 < f2
        alpha = knot(i-1);
        minf = f1;
        break;
    end
    % by the middle point for checking the active set
     x = x0 + 0.5 * (knot(i - 1) + knot(i))*p;
    Iu1 = x < 1e-15;
    x(Iu1) = 0;
    r = b - A * x;
    Iu2 = r < 1e-15;
    xb = x0;
    xb(Iu1) =0;
    rb = b - A * xb;
     % find the free components
     Ad=A(~Iu2,~Iu1)*p(~Iu1);
     % As the set Iu is full, NaN will appear
     if isempty(Ad)
         % The reason may be some error.
         alpha = knot(i - 1);
         minf = funmin(A,b,x0,p,alpha);
        return
     end
     % the derive Ap'(r-alpha*Ap)=0
     Ar=A(~Iu2,~Iu1)'*rb(~Iu2);
     % alpha = (Ap*r)^T/(Ap'*Ap)
     alpha = (p(~Iu1)'*Ar)/(Ad'*Ad);
     if alpha > knot(i)
         % it mean the interval is decreasing
         alpha = knot(i);
     else
         if alpha < 0
            alpha = knot(i-1);
         end
 %        minf = funmin(A,b,x0,p,alpha);
         break;
     end
end
minf = funmin(A,b,x0,p,alpha);

