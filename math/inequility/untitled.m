clear;clc
gamm = 0.6;
px1 = 0:0.2:1;
p1 = [px1' -px1'+1.5];fm1 = size(p1,1);
p2 = [px1' -px1'+0.5];fm2 = size(p2,1);
%plot(p1(:,1),p1(:,2),'r*',p2(:,1),p2(:,2),'+');
A1(1:fm1,:) = p1;
A1(fm1+1:fm1+fm2,:) = -p2;
b1 = [(1 + gamm)*ones(fm1,1);(gamm-1)*ones(fm2,1)];% 
x0=zeros(size(A1,2),1);maxIter = 10;
[xkh1,rkh,countFMh,countNWh,beginNWh,tfh,vkh,rkArrh]=han(x0,A1,b1,maxIter);
ww = [ xkh1';xkh1'; ];bb = [ gamm+1 ;1-gamm;];
d = lineData(ww , bb, [0,1], [0,1]);% ª≠∂‡Ãıœﬂ
%hold on 
%line(d(:,[1,2])',d(:,[3,4])')
[[xkh1' gamm+1]/norm(xkh1) norm(rkh);[xkh1' 1-gamm]/norm(xkh1) norm(A1'*rkh)]
