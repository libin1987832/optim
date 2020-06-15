clc
clear
addpath('../other')
addpath('../symsub')
addpath('../util')
addpath('../ours')
n=10000;
condtions=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6];
for i=1:1
    for j=1:1
        C=sprandsym(n,0.01,condtions(i),2);
        xs=sprandn(n,1,0.3);
        % xs must be nonnegative
        I1=(xs<0);
        xs(I1)=0;
        % count the number of active set
        I2=(xs>0);
        cs2=sum(I2);
        % generate q
        q=sprandn(n,1,0.3);
        qt=C*xs;
        q(xs>0)=-qt(xs>0);
        % nondegenency
        q(xs==0)=max(abs(qt))+0.1;
        
%         save('fpi1','C','xs','q','n');
        % subspace precise than our methods
        % load('fpi')
        
        x0=sparse(1:n,1:1,ones(n,1));
        nmax=10;
        
        max_iter = nmax;
        tol_rel  = 1e-5;
        tol_abs  = 1e-10;
        nf=5;
        
        
        % [xks,ress]=splitS(C,q,1.4,x0,10);
        disp(['mindig:' num2str(min(diag(C)))]);
        tic;[xk2,err,index2]=splitForlcp(x0,nmax,nf,C,q);f1=toc;
        tic;[xkpa,errpa,indexpa1,indexpa2]=PA(x0,nmax,nf,C,q);f2=toc;
        tic;[xkor,error,indexor,indexNor]=hybridorigin(x0,nmax,nf,C,q);f3=toc;
        tic;[xkpsor err iter flag convergence msg] = psor(C, q, x0, 1, max_iter*2, tol_rel, tol_abs, true);f4=toc;
        % || x-x^*||
        xspsor=norm(xs-xkpsor);
        xs2=norm(xs-xk2);
        xspa=norm(xs-xkpa);
        xsor=norm(xs-xkor);
        % ||KKT||
        [ress,fxs]=test_valid(C,q,xs);
        % [resss,fxss]=test_valid(C,q,xks);
        [res2,fx2]=test_valid(C,q,xkpsor);
        [res3,fx3]=test_valid(C,q,xk2);
        [res4,fx4]=test_valid(C,q,xkpa);
        [resor,fxor]=test_valid(C,q,xkor);
        disp([num2str(n),' ',num2str(condtions(i)),'&',num2str(res2),'&',num2str(max_iter),'&',...
            num2str(res3),'&',num2str(index2*nf+index2),'&',num2str(index2),'&',num2str(index2),'&',...
            num2str(resor),'&',num2str(indexor+indexNor),'&',num2str(indexor/nf),'&',num2str(indexNor)])
    end
end