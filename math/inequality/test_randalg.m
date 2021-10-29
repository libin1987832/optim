%% read data 
% [A,rows,cols,entries,rep,field,symm]=mmread('../util/well1033.mtx');
addpath('./dataInequality/');
addpath('./algorithmInequality/');
addpath('./randomized algorithm')

m = 1000;
n = 1000;
rangeMax = 2;
rangeMin = -2;

%    A = 2 * rand(m , n)-1;
%    b = 2 * rand(m , 1)-1;
%   % b=A*ones(n,1);
%    x0 = zeros(n , 1);
   A = -[1,-1;-1,-1;0,1];b=-[0;-1;0];x0=[-1;0];
   xkhe=[1/2;1/3];
   maxIter = 300;
   nf = 5;
   str = ['D','U','C','R','P'];
% for the solution so here

   % [xkh,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A,b,maxIter);
     [xkh,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A,b,x0,3,300,5,['D','HA']);
    [rk, rkh, dh, gh] = residual(A,b,xkh);
  %  fprintf('active:%d ,%g',vkh,dh);
    maxit = 50;
    xkh
    [xkacz,iterkacz,errorkacz,xAk,indexAk] = randomizedKaczmarzNE(A, b, x0,maxit,[],xkh);
    [xGS,iterGS,errorGS,xAg,indexAj] = randomizedGaussSeidelNE(A, b, x0,maxit,[],xkh);
    [xIn,iterIn,errorIn,xIng,indexAIn] = randomizedInexactNE(A, b, x0,maxit,[],xkh);
    
    
x = linspace(-0.5,1.5);
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
plot(xAk(1,:),xAk(2,:),'b+')
plot(xAg(1,:),xAg(2,:),'ro')
% plot(xIng(1,:),xIng(2,:),'g*')
line([0,1],[0,1]);
line([0,1],[1,0]);

indexAk
xAk
indexAj
xAg
indexAIn
xIng
 %%

beginp = 1;
figure
%maxIterA = 200;
h=semilogy(beginp:iterkacz,errorkacz,'b.');
h.LineStyle = '--';
hold on
h=semilogy(beginp:iterGS,errorGS);
h.LineStyle = '--';
legend('Kaczmarz','Gauss Seidel');
xlabel('the iterative numbers');
ylabel('the norm of the error');

%    end % for ration
%    end % for m