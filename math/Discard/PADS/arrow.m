function  arrow(P,V,color)
%��ά�ռ��л���ͷ
% ���룺P=[x0,y0],V=[a,b]
%����P��x0��y0��Ϊ��㣬�ԣ�x0+a��y0+b��Ϊ�յ㻭����ͷ
%���Խ�һ���޸�Ϊ��ά�ռ䵽��ͷ,��������PΪ��ʼ��,VΪ�յ�ļ�ͷͼ��

if nargin < 3 
    color = 'b'; 
end
x0 = P(1);y0 = P(2); 
a = V(1); b = V(2); 
l = max(norm(V), eps);
u = [x0 x0+a]; v = [y0 y0+b]; 
u1 = [x0 x0+a/2]; v1 = [y0 y0+b/2]; 
u2 = [x0+a/2 x0+a]; v2 = [y0+b/2 y0+b]; 
hchek = ishold;

%  plot(u,v,color)
 plot(u1,v1,color)
 hold on 
 plot(u2,v2,'-.')
%  plot(x0+a,y0+b,'rx')
% plot(x0+a,y0+b,'rx')
% h = l - min(.2*l, .2) ;v = min(.2*l/sqrt(3), .2/sqrt(3) );
% a1 = (a*h -b*v)/l*0.001;
% b1 = (b*h +a*v)/l*0.001; 
% 
% plot([x0+a1, x0+a], [y0+b1, y0+b], color) 
% 
% a2 = (a*h +b*v)/l; 
% b2 = (b*h -a*v)/l; 
% 
% plot([x0+a2, x0+a], [y0+b2, y0+b], color)
if hchek == 0 
    hold off 
end
