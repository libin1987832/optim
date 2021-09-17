%
% make_mtron
%
% Script to compile the mtron.c gateway routine. See readme.txt.
%


%
% On Mac OS X
%
% The following options were used on an Apple Powerbook, with OS X 10.4.7
% You may have to install g77 first. The make command for tron would be
%
%  make libs CC=gcc FC=g77 ARCH=darwin
%
mex mtron.c ./tron/src/tron/d-darwin.a ./tron/src/blas/d-darwin.a ./tron/src/icf/d-darwin.a ./tron/src/coloring/d-darwin.a -L'F:\msys64\usr\share\bash-completion\completions\gfortran'



%
% On Linux
% 
% The following options were used on some Fedora version. The make command
% for TRON was
% 
% make libs CC=gcc FC=g77 ARCH=linux
%
% mex mtron.c ./tron/src/tron/d-linux.a ./tron/src/blas/d-linux.a ./tron/src/icf/d-linux.a ./tron/src/coloring/d-linux.a -lg2c


%
% I have no idea how to compile TRON and mtron on Windows.
%
