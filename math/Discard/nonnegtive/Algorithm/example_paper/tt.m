A=[1,1/2;0,1/3];
x=[4;2];
xA=[];
for i=1:10
    xA=[xA x];
    y=A*x;
    x=y;
end
xA