% Matlab code RankDeficient.m
% For "Applied Numerical Linear Algebra", Figures 3.4 and 3.5
% Written by James Demmel, Feb 4, 1996
%                Modified, Jun 2, 1997
%
% Illustrate solution of Rank Deficient Least Squares Problem using
% Truncated SVD
%
% Given any m-by-n input matrix A and use supplied tolerance tol,
% we produce a matrix Al with a range of singular values from
% small to large. Then we compute a rank-deficient Ard by setting 
% the singular values of Al less than tol to zero. 
% Then we generate a rank-deficient least
% squares problem whose solution we know exactly as follows:
%
%  1) Compute the SVD Ard = U*S*V', where U and V have rnk columns,
%     where rnk is the rank of Ard.
%  2) Pick a random xrd of dimension rnk, let x=V*xrd, and 
%     let b = Ard*x (= U*S*xrd). Then x is the exactly, minimum norm
%     solution of the rank-deficient least squares problems 
%     min(norm(A*y-b,2)).
%  3) Add random noise dA of varying norm to A, and solve the perturbed
%     least squares problem using the truncated SVD, getting perturbed 
%     solution y.
%  4) Plot norm(y) and norm(y-x) versus norm(dA).
%
% Inputs:
%   dimensions m (number of rows) and n (number of columns), m >= n
%   type = 1 to computed matrix with floor(n/2) zero singular values;
%            for figure 3.4
%        = 2 to computed matrix with singular values roughtly geometrically
%            distributed between 1e-15 and 1; for figure 3.5
%   tol = tolerance, below which to set small singular value to 0
%
%
% Here are two ways to compute Al. One should be commented out.
m=20;
n=10;
type=1;
tol=1e-9;
A = randn(m,n);
if (type == 1)
%    Compute Al by zeroing out the last half of the columns of A:
    Al = randn(m,m) * ...
         [diag([ones(floor(n/2),1);zeros(ceil(n/2),1)]);zeros(m-n,n)] * ...
         randn(n,n);
end
if (type == 2)
%    Compute Al by multiplying columns by geometric sequence from 1e-16 to 1.
%    As a result Al will have singular values which form an approximate
%    geometric sequence in the same range.
     Al = A * diag(exp(log(10)*(0:n-1)*(-16/n)));
end
%
[Uorig,Sorig,Vorig]=svd(Al,0);
diag(Sorig);
rnk = length(find(diag(Sorig) >= tol));
disp(['rank(A) = ',int2str(rnk)])
U=Uorig(:,1:rnk);
S=Sorig(1:rnk,1:rnk);
V=Vorig(:,1:rnk);
Ard = U*S*V';
xrd = randn(rnk,1);
x = V * xrd;
b = Ard*x;
xnorm = norm(xrd);
nrm = sort([(10*ones(1,17)).^(-16:0),tol*(.1:.1:.9),tol*(.9:.01:1),...
             tol*(1:.1:1.9),tol*(2:10)]);
rnkdsav=[];
ynorm=[];
xmynorm=[];
for nn=nrm;
  dA = randn(m,n);
  dA = dA/norm(dA);
  dA = nn*dA;
  ArdpdA = Ard + dA;
  [Ud,Sd,Vd] = svd(ArdpdA);
  rnkd = length(find(diag(Sd) >= tol));
  rnkdsav=[rnkdsav;rnkd];
  Udd=Ud(:,1:rnkd);
  Sdd=Sd(1:rnkd,1:rnkd);
  Vdd=Vd(:,1:rnkd);
  y = Vdd*(Sdd\(Udd'*b));
  ynorm = [ynorm;norm(y)/xnorm];
  xmynorm = [xmynorm;norm(x-y)/xnorm];
end
figure(1),
hold off
%
subplot(2,1,1)
handl=semilogx(nrm',rnkdsav,'b'); set(handl,'LineWidth',2)
hold on
dS = min(max([diag(S)],1e-16),1);
handl=semilogx(dS,ones(size(dS)),'rx');set(handl,'MarkerSize',10)
axis([1e-16 1 0 min(m,n)]), grid
title(['Rank of perturbed A, original singular values are x''s, tol=',num2str(tol)])
xlabel('Norm of perturbation')
%
subplot(2,1,2)
handl=loglog(nrm',xmynorm,'b');set(handl,'LineWidth',2)
hold on
handl=semilogx(dS,ones(size(dS))*1e-15,'rx');set(handl,'MarkerSize',10)
axis([1e-16 1 1e-16 1]), grid
title('Norm(solution-perturbed solution)/norm(solution)')
xlabel('Norm of perturbation')