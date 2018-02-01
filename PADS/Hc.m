function [dy,dlamda]=Hc(H,c,x,l)
dy=H*x+c+l;
dy(find(dy>-1e5)&find(dy<1e5))=0;
dlamda=min(l,-1*x);