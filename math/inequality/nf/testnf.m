addpath('../dataInequality/');
            clc
            clear
for count=1:1
    for m=1000:1000:1000
        for ro=0.1:0.1:0.1

            % m = 2000;
            % n = 1800;
            rangeMax = 2;
            rangeMin = -2;
            n=floor(m*ro);     
            A = 2 * rand(m , n)-1;
            b = 2 * rand(m , 1)-1;
            x00 = zeros(n , 1);
            for k=3:10
                x0=x00;
                rk0 = b-A*x0;
                iter = 1500;
                [rk, rkN, nRkN, nGARk] = residual(A,b,x0);
                fprintf('residual %g gradient %g \n',nRkN,nGARk);
                [nf,iffind,xk,sumrkA]=generateNf(A,b,x0,iter,k);
                %[nf,iffind]=generateNf(A,b,x00,100);
                [rk, rkN, nRkN, nGARk] = residual(A,b,xk);
                fprintf('residual %g gradient %g ',nRkN,nGARk);
                sumrkA(1,iter-100:iter);
                iter2=100;
                name = [num2str(m) '_' num2str(n) '_' num2str(count) '_' num2str(k)];
                save([name,'.mat'],'A','b','x00','sumrkA');
                % plot(1:iter2,sumrkA(1:iter2));
                plot(1:iter2,abs(sumrkA(2:iter2+1)-sumrkA(1:iter2)));
                saveas(gcf,['D:\matlab',name,'.jpg']);
            end
        end
    end
end
% rkA=zeros(m,iter);
% sumrkA=zeros(1,iter);
% rkA(:,1)=rk0;
% sumrkA=sum(rk0>0);
% for i = 2:iter
% [xk,rk]=Lei(x0,A,b,3,rk0);
% x0=xk;
% rk0=rk;
% rkA(:,i)=rk;
% sumrkA(:,i)=sum(rk>0);
% end
% sumrkA