function [rkD,dD,gD,sumiter]=printLatex(xkD,A,b,itersm,nf,tfD,type)
        rkD=b-A*xkD;
        rkD(rkD<0)=0;
        dD=norm(rkD);
        gD=norm(A'*rkD);
%        beginN=find(itersm>0);
        sumiter=sum(itersm>0);
%         fprintA(2*i-1:2*i)=[gD,tfD];
%         if isempty(beginN)
%             beginN=0;
%         end
      %  record(i,:)=[dD,gD,iter*nf,sumiter,tfD];
% print for tex      
        % fprintf('%s $ %d \\times %d $ & %g & %g & %g & %g & %g & %g\n',type,m,n,dD,gD,tfD,iter*nf,sumiter,beginN(1));
       fprintf('& %s & %g & %g & (%d,%d)  & %g \\\\\n',[type,'HA'],dD,gD,iter*nf,sumiter,tfD);     
  
