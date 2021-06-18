clear 
clc
gg=semilogx((1e-1).^(1:10),1:10,'b*-',(1e-1).^(1:10),1:10,'g*-')
% set(gca,'XTickLabel',(1e-1).^(10:1));
set(gca,'XDir','reverse')