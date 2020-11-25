clc
clear
addpath('../FM')
addpath('../IFM')
addpath('../strategies')
addpath('../hybrid')
addpath('../linearsolve')
Arecord=[];
maxIter=500;
nmax=500;


m=1000;

for m=[6]
    n=3;
    batchArr=zeros(1,3*10);
    eigenA=zeros(2,10);
    ArrA=cell(10,1);
    for batch=1:10
         A=2*rand(m,n)-1;
       %A = 2*randn(m,m) * ...
        % [diag([ones(floor(n/2),1)*100;ones(ceil(n/2),1)])*10;ones(m-n,n)*10] * ...
        % randn(n,n)-1;
%         eA=eig(A'*A);
        b=2*rand(m,1)-1;
        x0=zeros(n,1);
        xr=rand(n,1);
        xs=-1;
        
        [xkhan,rkhan,countFMhan,countNWhan,beginNWhan,tfhan,vkhan,xkAhan]=han(x0,A,b,maxIter); 
        dhan=norm(xkhan-xs);
        xs=xkhan;
         rkhan=b-A*xkhan;
         rkhan(rkhan<0)=0;
         ghan=norm(A'*rkhan);
         fprintf('han$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,dhan,ghan,tfhan,countFMhan,countNWhan);
        
        [xkR,xkR2,countFR,countNWR,bNWR,tfR,vkR]=residualR(x0,A,b,maxIter);
        dR=norm(xkR-xs);
        rkR=b-A*xkR;
        rkR(rkR<0)=0;
        gR=norm(A'*rkR);
        dR2=norm(xkR2-xs);
        rkR2=b-A*xkR2;
        rkR2(rkR2<0)=0;
        gR2=norm(A'*rkR2);
        fprintf('resdual$ %d \\times %d $ & %g & %g & %4.2f & %d & %d &\n',m,n,gR2,gR,tfR,countFR,countNWR);
        
        [xkA,rkA,countFA,countNA,bNWA,tfA,vkA,Arr,rkrr]=als(x0,A,b,maxIter);
        dA=norm(xkA-xs);
        rkA=b-A*xkA;
        rkA(rkA<0)=0;
        gA=norm(A'*rkA);
        if ~isempty(Arr)
        fprintf('ALS$ %d \\times %d $ & %g & %g & %4.2f & %d & %d & %d\n',m,n,dA,gA,tfA,countFA,countNA,Arr(1,end));
        batchArr(1,3*batch-2:3*batch)=[gA,tfA,Arr(1,end)];
        %eigenA(:,batch)=[max(eA);min(eA)];
        ArrA{batch}=Arr;
        end
        %          
    end
    Arecord=[Arecord;m n batchArr];
end
% save(['ff' num2str(batch) '.mat'],'Arecord')
%%
nf=Arecord(:,5:3:32);
time=Arecord(:,4:3:32);
accuracy=Arecord(:,3:3:32);

%%
nfms=[nf mean(nf,2) std(nf,0,2)];
timems=[time mean(time,2) std(time,0,2)];
accuracyms=[accuracy mean(accuracy,2) std(accuracy,0,2)];