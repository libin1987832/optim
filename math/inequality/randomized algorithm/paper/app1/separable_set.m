clear
clc
gamm = 0.6;
n = 200000;
w = [1,1]+rand(1,2);
b = 1+rand(1);
x = rand( n , 2 );
%% triangle
w1 = [1,1];w2 = [1,-1];
b1 = 0.2;b2 = 1;
y1 = (x*w1'-(b1+1)*ones(n,1) < 0 & x*w2'-b1*ones(n,1) < 0);
y2 = (x*w1'-(1-b1)*ones(n,1) > 0 & x*w2'+b1*ones(n,1) > 0);
p3 = x(y1&y2,:);
row = size(p3,1)/2;
p1 = x(y1&~y2,:);
p2 = x(~y1&y2,:);
p1 = [p1;p3(1:row,:)];
p2 = [p2;p3(row+1:end,:)];


%%
fm1 = size(p1,1);
fm2 = size(p2,1);
A1(1:fm1,:) = p1;
A1(fm1+1:fm1+fm2,:) = -p2;
A2 = [[p1 -1*ones(fm1,1)];[-p2 ones(fm2,1)]];
% b1 = [(1 + gamm)*ones(fm1,1);(-1-gamm)*ones(fm2,1)];% 
b2 = [ones(fm1,1);ones(fm2,1)];
source_train = [A1(1:fm1,:);-A1(fm1+1:fm1+fm2,:)];
label_train = [ones(fm1,1);-1*ones(fm2,1)];
%label_train = round(rand(fm1+fm2,1));


x0=zeros(size(A1,2),1);
maxIterA = 0;
maxIter = 2000;
nf = 5;
str = ['D','C','R','P','H'];
steplength = 0;%1/(max(eig(ATA))+0.0001);
% A=A2;
% 
tol=-1;
 x_exact = [] ;
 debug = 0;
% x1=[x0;0];
% t=clock;
% [x2,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A2,b2,[x0;0],maxIter,nf,'RHA');
% xst(2)=etime(clock,t);
% t=clock;
% % % [x3,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A,b,x1,maxIter,nf,'PHA');
% % xst(3)=etime(clock,t);
% % t=clock;
% % % [x4,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=hybridA(A,b,x1,maxIter,nf,'CHA');
% % xst(4)=etime(clock,t);
% % 
%  xs=[x2 x2 x2 x2];
% % 
% str=['P','R','P','C'];
% for i = 2:2
% r = b - A * xs(:,i);
% r(r<0) = 0;
% r_GS = norm(r);
% g_GS = norm(A'*r);
% fprintf('& %sHA  &%g & %g & %g \\\\\n', str(i),r_GS, g_GS, xst(i));
% end





%% GuassSeidel
maxit_Rand =20000;
t=clock;
 [x_GS,iter_GS,error_GS,xA_GS,index_GS] = GuassSeidelNE(A2,b2,[x0;0],2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b2 - A2 * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A2'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'GuassSeidel', r_GS, g_GS, iter_GS, tf_GS);


%% simpleGuassSeidel
% maxit_Rand =1000000;
t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = simpleGuassSeidelNE(A2,b2,[x0;0],2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b2 - A2 * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A2'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'simpleGuassSeidel', r_GS, g_GS, iter_GS, tf_GS);

%% randGuassSeidel
% maxit_Rand =630000;
t=clock;
[x_GS,iter_GS,error_GS,xA_GS,index_GS] = randGuassSeidelNE(A2,b2,[x0;0],2.0,maxit_Rand,tol,x_exact,debug);
tf_GS=etime(clock,t);
r = b2 - A2 * x_GS;
r(r<0) = 0;
r_GS = norm(r);
g_GS = norm(A2'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n', 'randGuassSeidel', r_GS, g_GS, iter_GS, tf_GS);



%% 参数的设定
 maxit_IFM = 2000;
maxit_LSQR =3;

t=clock;
[x_IFM,iter_IFM,error_IFM,xA_IFM,index_IFM] = IFM(A2,b2,[x0;0], maxit_IFM, maxit_LSQR ,tol, x_exact,debug);
tf_IFM=etime(clock,t);
r = b2 - A2 * x_IFM;
r(r<0) = 0;
r_IFM = norm(r);
g_IFM = norm(A2'*r);
fprintf('& %s & %g & %g & %d & %g \\\\\n','IFM', r_IFM, g_IFM,iter_IFM,tf_IFM);

% FM
% maxit =100;
% 
% alpha=1;
% maxit_gs=n;
% t=clock;
% [x_FM,iter_FM,error_k,iter_dFM,index_k] = DFM(A2,b2,[x0;0], maxit_IFM,alpha,maxit_gs,tol, x_exact,debug);
% tf_FM=etime(clock,t);
% r = b - A * x_FM;
% r(r<0) = 0;
% r_FM = norm(r);
% g_FM = norm(A'*r);
% fprintf('& %s & %g & %g & %d & %g \\\\\n','FM', r_FM, g_FM,iter_FM,tf_FM);






%  for i=1:5
%         type = str(i);   
%         if i<5
%         [xkD,flag,relres,iter,resvec,arvec,itersm,tfD]=hybridA(A2,b2,[x0;0],steplength,maxIter,nf,[type,'HA']);
%         resvec = arvec;
%         else
%         [xkD,rkh,countFMh,~,itersm,tfD,vkh,rkArrh]=han([x0;0],A2,b2,maxIter); 
%         end%fprintf('& %s & %g \\\\\n',[type,'HA'],tfD);
%         % plot picture
%         A=A2;
%         b=b2;
%         resvec = arvec;
% % check the solution
%         rkD=b-A*xkD;
%         rkD(rkD<0)=0;
%         dD=norm(rkD);
%         gD=norm(A'*rkD);
%         beginN=find(itersm>0);
%         sumiter=sum(itersm>0);
%         fprintA(2*i-1:2*i)=[gD,tfD];
%         if isempty(beginN)
%             beginN=0;
%         end
%         record(i,:)=[dD,gD,iter*nf,sumiter,tfD];
% % print for tex      
%         % fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
%        fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);
%  end
% [xkh2,rkh,countFMh,countNWh,itersm,tfD,vkh,rkArrh]=han([x0;0],A2,b2,maxIter); 


% ww = [w; wsvn'; xkh2([1,2])';xkh2([1,2])'; xkh2([1,2])' ];
% bb = [b; -bsvn; xkh2(3)+1 ;xkh2(3)-1;xkh2(3)];
% display = [3,4,5];

% 画多条线
% d = lineData(ww(display,:) , bb(display,:) , [0,1], [0,1]);
% plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');
% hold on 
% line(d(:,[1,2])',d(:,[3,4])')
% [[w b]/norm(w);[wsvn' -bsvn]/norm(wsvn);[xkh2([1,2])' xkh2(3)]/norm(xkh2([1,2]))];

 
 
 
  %% plot picture
% type=['r','c','k','g'];
% typet=['+','o','v','s'];
% beginp = 1;
% figure
% %maxIterA = 70;
% h=semilogy(beginp:maxIterA,iterA(1,beginp:maxIterA),'bx');
% h.LineStyle = '--';
% hold on
% for i=2:4
%     h=semilogy(beginp:maxIterA,iterA(i,beginp:maxIterA),[type(i) typet(i)]);
%     h.LineStyle = '--';
% end
% legend('DHA','GHA','RHA','PHA');
% xlabel('Iteration Number');
% ylabel('the norm of the gradient');
