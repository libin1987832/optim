x1 = imread('untitled.tif');
x2 = imread('untitled2.tif');
x3 = imread('untitled3.tif');
x=[x1 x2 x3];
imshow(x)
% h1=openfig('untitled.fig','reuse');
% h2=openfig('untitled2.fig','reuse');
% h3=openfig('untitled3.fig','reuse');
% figure
% s1=subplot(1,3,1);
% copyobj(get(get(h1,'Children'),'Children'),s1)
% s2=subplot(1,3,2);
% copyobj(get(get(h2,'Children'),'Children'),s2)
% s3=subplot(1,3,3);
% copyobj(get(get(h3,'Children'),'Children'),s3)