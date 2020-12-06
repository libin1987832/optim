function [A,b,x0] = readData(type)
switch type
    case 1
        m1 = 5;
        m2 = 5; n = 7;
        A1=abs(sprand(m1,n,0.1,1/100))/100;
        A2=abs(sprand(m2,n,0.1,1/100))/100;
        %A1=rand(m1,n) + 1 ;
        %A2=rand(m2,n) + 1;
        b1=rand(m1,1);
        b2=rand(m2,1);
        A=[A1;-A2];
        b=[b1;-b2];
        x0=sparse(ones(n,1));
    case 2
         load('test')
         m1 = 5;m2=5;n=7;
    otherwise
end