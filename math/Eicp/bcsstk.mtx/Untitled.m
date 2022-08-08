clear;
clc;
filename='bcsstk02.mtx';
[A,rows,cols,entries,rep,field,symm]=mmread(filename);

eig(A)
