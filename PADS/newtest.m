%produce H c 
%gx globe optimization lx local optimization n dimesion I index set realA
%active set sep sept hydral plane allSet all information 

% H=[4 5 -5;5 9 -5;-5 -5 7];
% c=[2;1;-3];
% lx=[-0.5;0;0];
% n=3;

function [allSet,cycle,realcycle,realA,sep]=newtest(H,c,lx)
n=size(lx,1);
gx=-1*inv(H)*c;
I=1:n;
realA=I(lx==0);
sep=I(gx>0);
nsep=setdiff(I,sep);
%disp(['globf:',num2str(-0.5*c'*inv(H)*c),' real active:',num2str(realA),' sep:',num2str(sep)]);
allSet=cell(2^n,10);
maxX=-inf;minX=inf;maxT=-inf;minT=inf;
for i=0:2^n-1
    bb=bitget(uint8(i),1:1:n);
    A=I(bb>0);
    [x,z,Ak,nextA,f,kktz,kkta]=PADSA(H,c,A,n);
    %x next iteration x .z dual variable,A input work set,Ak next work set 
    %f object function kkt KKT norm i current work set number(the binary number is work set index)
    %next next work set number
    allSet(i+1,:)={x,z,A,Ak,f,kktz+kkta,i,nextA,kktz,kkta};
%     disp(['active set:', num2str(A) ,' next active set:' , num2str(Ak) , ...
%         ' fun:' , num2str(f) , ' kkt:' , num2str(kktz+kkta),...
%         ' kktz:',num2str(kktz),' kkta:',num2str(kkta)]);
     disp(['kktz:',num2str(kktz),'       kkta:',num2str(kkta)]);
end
%every step record link,first is number vector,second is whether is cycle,third
%is string link. realcycle the same cycle only record once.
cycle=cell(2^n,3);
realcycle=cell(1);
for i=0:2^n-1
    cycle{i+1,1}=[i];
    cycle{i+1,3}=['[',num2str(allSet{i+1,3}),']-->'];
    fprintf('%s','active link:[');
    diffActive(allSet{i+1,3},realA,sep);
    fprintf('%s',']-->');
    nexti=i;
    %cur find link until cycle occur.
    while 1
        nexti=allSet{nexti+1,8};
        index=ismember(cycle{i+1,1},nexti);
        cycle{i+1,1}=[cycle{i+1,1},nexti];
        cycle{i+1,3}=[cycle{i+1,3},'[',num2str(allSet{nexti+1,3}),']-->'];
        fprintf('%s','[');
        diffActive(allSet{nexti+1,3},realA,sep);
        fprintf('%s',']');
        sindex=size(index,2);
        %find a cycle,finish fill cycle if exits and the cycle if new 
        if find(index>0) > 0
             if index(1,sindex)
                cycle{i+1,2}='Y';
                cycle{i+1,3}=[cycle{i+1,3},'YLOCLA'];
                fprintf('                             %s\n','YLOCLA');
             else
                cycle{i+1,2}='N';
                cycle{i+1,3}=[cycle{i+1,3},'NLOCLA'];
                fprintf('\n');
             end
            cv=cycle{i+1,1};
            tRealcy=cv(1,find(index==1):1:sindex);
            %sort only compared,if there is exit.
            tRealcs=sort(tRealcy);
            stTtRealcs=size(tRealcs,2);
            stRealcs=size(realcycle,1);
            for j=1:stRealcs
                % the cycle is exits so break.
                if size(realcycle{j,1},2)==stTtRealcs && sum(realcycle{j,1}-tRealcs) < 1
                    break;
                end
                if j==stRealcs
                 realcycle{stRealcs+1,1}=tRealcs;
                 realcycle{stRealcs+1,2}=[];
                 realcycle{stRealcs+1,3}=tRealcy;
                 for k=1:size(tRealcs,2)-1
                     realcycle{stRealcs+1,2}=[realcycle{stRealcs+1,2},'[',...
                         num2str(allSet{tRealcy(k)+1,3}),']-->'];
                 end
                 realcycle{stRealcs+1,2}=[realcycle{stRealcs+1,2},'[',...
                     num2str(allSet{tRealcy(stTtRealcs)+1,3}),']'];
                end
            end
            break;
        end
        fprintf('%s','-->');
    end
%     tSet=intersect(allSet{i+1,3},realA);
%     SSet=intersect(allSet{i+1,3},sep);
    %disp(['active link:',cycle{i+1,3}]);
%    cprintf('text','active link:%s',cycle{i+1,3});
%     cprintf('key',' real active set:%s ',num2str(tSet));
%     cprintf('cyan',' sep set:%s\n',num2str(SSet));
end
    disp(['globf:',num2str(-0.5*c'*inv(H)*c),' real active(A):',num2str(realA),' sep(S):',num2str(sep)]);
for i=2:size(realcycle,1)
    disp(['cycle link:',realcycle{i,2}]);
end
% maxf=max([allSet{:,5}])+0.0001;
% minf=min([allSet{:,5}])-0.0001;
% maxk=max([allSet{:,6}])+0.0001;
% mink=min([allSet{:,6}])-0.0001;
% axis([min([allSet{:,6}]),max([allSet{:,6}]),min([allSet{:,5}]),max([allSet{:,5}])]);
% for i=1:2^n
%    plot(allSet{i,5},allSet{i,6},'*');
%    text(allSet{i,5},allSet{i,6},strcat('[ ',strcat(num2str(allSet{i,3}),' ]')));
%     arrowX=[allSet{i,5}-minf,allSet{allSet{i,8}+1,5}-minf]./(maxf-minf);
%     arrowY=[allSet{i,6}-mink,allSet{allSet{i,8}+1,6}-mink]./(maxk-mink);
%     %ah=annotation('arrow',arrowX,arrowY,'Color','r');
%     arrow([allSet{i,5},allSet{i,6}],[allSet{allSet{i,8}+1,5},allSet{allSet{i,8}+1,6}]-[allSet{i,5},allSet{i,6}]);
%    hold on
% end
% grid on;
