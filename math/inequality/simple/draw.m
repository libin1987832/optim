 figure
 x1=imread('FM.tif');
 % subplot(221);
 % imshow(x);
 x2=imread('IFM.tif');
 
 % subplot(222);
 % imshow(x);
 x3=imread('HAN.tif');
 % subplot(223);
 % imshow(x);
x4=imread('DHA2.tif');

 % subplot(224);
 % imshow(x);
 
 x5 = [x1 x2 ;
     x3 x4];
%  imshow(x5)
%   imwrite(x5,'tt.eps')
  imshow(imresize(x5,0.9))
 % copyobj(get(get(h1,'Children'),'Children'),s1)
 % s2=subplot(212);
 % copyobj(get(get(h2,'Children'),'Children'),s2)
 % s3=subplot(2,2,3);
 % copyobj(get(get(h3,'Children'),'Children'),s3)