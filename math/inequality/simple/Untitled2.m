


% h1=openfig('FM.fig','reuse');
% h2=openfig('IFM.fig','reuse');
% h3=openfig('HAN.fig','reuse');
% h4=openfig('DHA.fig','reuse');
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
x4=imread('DAX2.tif');
% subplot(224);
% imshow(x);

x5 = [x1 x2 ;
    x3 x4];
imshow(x5)
% copyobj(get(get(h1,'Children'),'Children'),s1)
% s2=subplot(212);
% copyobj(get(get(h2,'Children'),'Children'),s2)
% s3=subplot(2,2,3);
% copyobj(get(get(h3,'Children'),'Children'),s3)
% s4=subplot(2,2,4);
% copyobj(get(get(h4,'Children'),'Children'),s4)

% h1=openfig('FM.fig','reuse');
% h2=openfig('IFM.fig','reuse');
% h3=openfig('HAN.fig','reuse');
% h4=openfig('DHA.fig','reuse');
% figure
% s1=subplot(211);
% copyobj(get(get(h1,'Children'),'Children'),s1)
% s2=subplot(212);
% copyobj(get(get(h2,'Children'),'Children'),s2)
% s3=subplot(2,2,3);
% copyobj(get(get(h3,'Children'),'Children'),s3)
% s4=subplot(2,2,4);
% copyobj(get(get(h4,'Children'),'Children'),s4)
% h1=openfig('FM.fig','reuse');
% h2=openfig('IFM.fig','reuse');
% h3=openfig('HAN.fig','reuse');
% h4=openfig('DHA.fig','reuse');
% figure
% s1=subplot(211);
% copyobj(get(get(h1,'Children'),'Children'),s1)
% s2=subplot(212);
% copyobj(get(get(h2,'Children'),'Children'),s2)
% s3=subplot(2,2,3);
% copyobj(get(get(h3,'Children'),'Children'),s3)
% s4=subplot(2,2,4);
% copyobj(get(get(h4,'Children'),'Children'),s4)

% h1=openfig('FM.fig','reuse');
% h2=openfig('IFM.fig','reuse');
% figure
% s1=subplot(211);
% cc1=get(get(h1,'Children'),'Children');
% copyobj(cc1{2},s1)
% s2=subplot(212);
% cc2=get(get(h2,'Children'),'Children');
% copyobj(cc2,s2)