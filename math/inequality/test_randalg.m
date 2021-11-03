%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
addpath('./dataInequality/');
addpath('./algorithmInequality/');
addpath('./randomized algorithm')

m = 10000;
n = 100;
rangeMax = 2;
rangeMin = -2;

%    A = 2 * rand(m , n)-1;
%    b = 2 * rand(m , 1)-1;
%   % b=A*ones(n,1);
%    x0 = zeros(n , 1);
    A = [1,-1;-1,-1;0,1];b=[0;-1;0];x0=[-1;0];
    xkhe=[1/2;1/3];
   maxIter = 20;
   nf = 5;
   str = ['D','U','C','R','P'];
% for the solution so here

   % [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
    %[xkh,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,3,300,5,['R','HA']);
   % [rk, rkh, dh, gh] = residual(A,b,xkh);
  %  fprintf('active:%d ,%g',vkh,dh);
    maxit = 10;
    xkh;
    norm(gh)
    mutiple = 10;
   [xIn,iterIn,errorIn,xIng,indexAIn] = randomizedInexactNE(A, b, x0,maxit*mutiple,[],xkh);
     [xkacz,iterkacz,errorkacz,xAk,indexAk] = randomizedKaczmarzNE(A, b, x0,maxit,[],xkh);
   % [xkacz,iterkacz,errorkacz,xAk,indexAk] = IFM(A, b, x0,maxit,[],xkh);
    [xGS,iterGS,errorGS,xAg,indexAj] = randomizedGaussSeidelNE(A, b, x0,maxit*mutiple,[],xkh);

    
    
x = linspace(-1,1.5);
y = linspace(0,1);
[X,Y] = meshgrid(x,y);
XZ = repmat(X,1,1,3);
YZ = repmat(Y,1,1,3);
ba1 = reshape(A(:,1),1,1,3);
ba2 = reshape(A(:,2),1,1,3);
br = reshape(b,1,1,3);
z = bsxfun(@times,ba1,XZ )+bsxfun(@times,ba2,YZ );
z = bsxfun(@minus,br,z);
z(z<0)=0;
z=0.5*z.^2;
Z=squeeze(sum(z,3));
figure
contour(X,Y,Z)
hold on
plot(xAk(1,4),xAk(2,4),'b+')
% plot(xAg(1,:),xAg(2,:),'ro')
% plot(xIng(1,:),xIng(2,:),'g*')
line([0,1],[0,1]);
line([0,1],[1,0]);

% indexAk
% xAk
% indexAj
% xAg
% indexAIn
% xIng
 %%

% beginp = 1;
% figure
% maxIterA = 200;
% count = 50;
% jump=floor(iterkacz/count)+1;
% h=semilogy((beginp:jump:iterkacz)*mutiple,errorkacz(beginp:jump:iterkacz),'g.');
% h.LineStyle = '--';
% hold on
% jump=floor(iterGS/count);
% h=semilogy(beginp:jump:iterGS,errorGS(beginp:jump:iterGS),'r+');
% h.LineStyle = '--';
% jump=floor(iterIn/count);
% h=semilogy(beginp:jump:iterIn,errorIn(beginp:jump:iterIn),'b*');
% h.LineStyle = '--';
% legend('IFM','Gauss Seidel','Inexact');
% % legend('Gauss Seidel','Inexact');
% xlabel('the iterative numbers');
% ylabel('the norm of the error');

%    end % for ration
%    end % for m