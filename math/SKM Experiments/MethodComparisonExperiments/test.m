diary('output.txt')
diary on
time = [];
fun = @(x) 0;
A = randn(1000,100);
b = A*ones(100,1) + abs(randn(1000,1));

tstart = cputime;
x = fmincon(fun,zeros(100,1),A,b);
time(1,1) = cputime - tstart

tstart = cputime;
x = SKM(A,b,zeros(100,1),1,10^(-6),10000,50);
time(1,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set');
tstart = cputime;
x = fmincon(fun,zeros(100,1),A,b,[],[],[],[],[],options);
time(1,3) = cputime - tstart

load('lp_adlittle.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
% column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*10*ones(n,1) - d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',10^(-2),'TolX',0,'TolFun',10^(-2),'MaxFunEvals',1000,'Display','iter');
tstart = cputime;
x = fmincon(fun,10*ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(2,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,10*ones(n,1),1.2,10^(-2)*initres,10000,30);
time(2,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',10^(-2),'TolFun',10^(-2),'Display','iter');
tstart = cputime;
x = fmincon(fun,10*ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(2,3) = cputime - tstart

load('lp_agg.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);

options = optimoptions('fmincon','Algorithm','interior-point','TolCon',10^(-2),'TolFun',10^(-2),'MaxFunEvals',10000,'Display','iter');%10^(-2)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(3,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1,10^(-2)*initres,10000,100);
time(3,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',10^(-2),'TolFun',10^(-2),'Display','iter');
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(3,3) = cputime - tstart

load('lp_bandm.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n);
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',10^(-2),'TolFun',10^(-2),'MaxFunEvals',100000);%10^(-2)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(4,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1.2,10^(-2)*initres,10000,100);
time(4,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',10^(-2),'TolFun',10^(-2),'Display','iter');
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(4,3) = cputime - tstart

load('lp_blend.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',10^(-3),'TolFun',10^(-3),'MaxFunEvals',100000);%10^(-3)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(5,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1.6,10^(-3)*initres,10000,250);
time(5,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',10^(-3),'TolFun',10^(-3));
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(5,3) = cputime - tstart

load('lp_brandy.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',5*10^(-2),'TolFun',5*10^(-2),'MaxFunEvals',100000);%5*10^(-2)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(6,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1,5*10^(-2)*initres,10000,20);
time(6,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',5*10^(-2),'TolFun',5*10^(-2));
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(6,3) = cputime - tstart

load('lp_degen2.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',10^(-2),'TolFun',10^(-2),'MaxFunEvals',100000);%10^(-2)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(7,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1.4,10^(-2)*initres,10000,100);
time(7,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',10^(-2),'TolFun',10^(-2),'Display','iter');
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(7,3) = cputime - tstart

load('lp_finnis.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',5*10^(-2),'TolFun',5*10^(-2),'MaxFunEvals',100000,'Display','iter');%5*10^(-2)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(8,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1,5*10^(-2)*initres,10000,50);
time(8,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',5*10^(-2),'TolFun',5*10^(-2),'Display','iter');
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(8,3) = cputime - tstart

load('lp_recipe.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',2*10^(-3),'TolFun',2*10^(-3),'MaxFunEvals',100000);%2*10^(-3)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(9,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1.2,2*10^(-3)*initres,10000,30);
time(9,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',2*10^(-3),'TolFun',2*10^(-3));
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(9,3) = cputime - tstart

load('lp_scorpion.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',5*10^(-3),'TolFun',5*10^(-3),'MaxFunEvals',100000);%5*10^(-3)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(10,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1.6,5*10^(-3)*initres,10000,200);
time(10,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',5*10^(-3),'TolFun',5*10^(-3));
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(10,3) = cputime - tstart

load('lp_stocfor1.mat');
A = full(Problem.A);
b = Problem.b;
fun = @(x) (Problem.aux.c)'*x;
m = size(C,1);
n = size(C,2);
%column vector with norm of each row
norm_rowsA=sqrt(sum((C).^2,2));

%normalize vector b
d(1:m) = d(1:m)./norm_rowsA;   

%normalize matrix A
C(1:m,:) = C(1:m,:)./repmat(norm_rowsA,1,n); 
C = full(C);
initres = max(C*ones(n,1)-d);
options = optimoptions('fmincon','Algorithm','interior-point','TolCon',10^(-1),'TolFun',10^(-1),'MaxFunEvals',100000);%10^(-1)
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(11,1) = cputime - tstart

tstart = cputime;
x = SKM(C,d,ones(n,1),1.4,10^(-1)*initres,10000,50);
time(11,2) = cputime - tstart

options = optimoptions('fmincon','Algorithm','active-set','TolCon',10^(-1),'TolFun',10^(-1));
tstart = cputime;
x = fmincon(fun,ones(n,1),[],[],A,b,Problem.aux.lo,Problem.aux.hi,[],options);
time(11,3) = cputime - tstart

diary off