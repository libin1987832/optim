function [I,S,L,E]=prerepeat(A,c.y0)
S=c-A'*y0;
I=(S>=0)
[E,L]=eig(A*A)'