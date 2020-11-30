% interv [0,1] for choose a 
function [alpha, minf, knot,retcode] = piecewiseNnls(A,b,x0,p,tol,maxits)
if nargin<5, tol = 1e-15; maxits = 1e4; end
[m,n] = size(A);
nonzerou = p < 0;
arange = -x0(nonzerou)./p(nonzerou);
aranges = sort(arange(arange > tol & arange < 1));
aranges = [0;aranges;1];
aranges = unique(aranges);
all = zeros(m*n,1);
[] = arrayfun(@(x) xa = x0 + x.a; xa(xa<tol)=0; r =b-A*xa;r)
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
    ai=r./ap;
    alength = aranges( i ) - aranges( i -1);
    as=sort(ai(ai > 0 & alength >ai));
    ass = aranges(i-1) + as;
    ass = [aranges( i-1 ) ass' aranges(i)]; 
    all = [all,ass];
end
knot = unique(all);
for i = 2 : length(knot)
    f1 = funmin(A,b,x0,p,knot( i - 1 ));
    f2 = funmin(A,b,x0,p,knot( i ));
    if f1 < f2
        alpha = knot(i-1);
        minf = f1;
        break;
    end
     x = x0 + 0.5 * (knot(i - 1) + knot(i))*p;
    Iu1 = x < 1e-15;
    x(Iu1) = 0;
    r = b - A * x;
    Iu2 = r < 1e-15;
    xb = x0;
    xb(Iu1) =0;
    rb = b - A * xb;
      % Ap
     Ad=A(~Iu2,~Iu1)*p(~Iu1);
     if isempty(Ad)
         alpha = knot(i - 1);
         minf = funmin(A,b,x0,p,alpha);
        return
     end
     % the derive Ap'(r-alpha*Ap)=0
     Ar=A(~Iu2,~Iu1)'*rb(~Iu2);
     % alpha = (Ap*r)^T/(Ap'*Ap)
     alpha = (p(~Iu1)'*Ar)/(Ad'*Ad);
     if alpha > knot(i)
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

