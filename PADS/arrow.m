function  arrow(P,V,color)
%二维空间中画箭头
% 输入：P=[x0,y0],V=[a,b]
%将以P（x0，y0）为起点，以（x0+a，y0+b）为终点画出箭头
%可以进一步修改为三维空间到箭头,或者是以P为起始点,V为终点的箭头图像

if nargin < 3 
    color = 'b'; 
end
x0 = P(1);y0 = P(2); 
a = V(1); b = V(2); 
l = max(norm(V), eps);
u = [x0 x0+a]; v = [y0 y0+b]; 
hchek = ishold;

plot(u,v,color)
hold on 
plot(x0+a,y0+b,'rx')
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
