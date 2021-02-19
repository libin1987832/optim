classdef parameter
    % precision
    properties
        % hybrid
        tol
        % simple
        % sm
        lsqrTol
        % strategy
    end
    % parameter
    properties
        % hybrid
        maxIter
        nf
        type
        % simple
        steplengthOrk
        % sm
        lsqrIter
        smIsqmIter
        % strategy
        % dax
        % contract
        con1
        con2
        % reduction
        diff
        % prediction
        eIter
        eIter_num
    end
    methods
%picture_paper 1000 800 save('pictureconsistent.mat','A','b','x0')
        %  \hline 
%  \multirow{4}{*}{$ 1000\times 800 $}& DHA & 9.04455e-16 & 1.39461e-14 & (130,0,0)  & 0.906 \\
% & CHA & 1.23021e-15 & 1.8318e-14 & (20,2,0)  & 0.221 \\
% & RHA & 2.61829e-15 & 3.77841e-14 & (25,1,1)  & 0.227 \\
% & PHA & 3.11251e-15 & 4.48269e-14 & (15,1,1)  & 0.173 \\       
  function param = parameter()
            param.tol = 1e-15;
            param.maxIter = 200;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 3;
            param.con1 = 0.95;
            param.con2 =0.6;
            param.diff = 100*eps;%consistent 10*eps
            param.eIter = 5;
            param.eIter_num = 0.8;
        end

        
        
        function param = parameter_()
            param.tol = 1e-15;
            param.maxIter = 600;
            param.nf = 3;
            param.steplengthOrk = 10;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.3;
            param.diff = 1e-10;
            param.eIter = 5;
            param.eIter_num = 0.95;
        end

%data  [A,rows,cols,entries,rep,field,symm]=mmread('well1850.mtx');
        %         \multirow{4}{*}{$ 1850\times 712 $}& IFM & 9.59166 & 4.94528e-16 & (1500,0)  & 0.196 \\
% & DHA & 9.59166 & 9.7861e-16 & (516,0)  & 0.174 \\
% & CHA & 9.59166 & 1.16511e-13 & (9,2)  & 0.111 \\
% & RHA & 9.59166 & 1.96507e-13 & (57,6)  & 0.117 \\
% & PHA & 9.59166 & 4.05799e-10 & (63,9)  & 0.148 \\
          function param = parameter_well1850()
            param.tol = 1e-15;
            param.maxIter = 600;
            param.nf = 3;
            param.steplengthOrk = 10;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.3;
            param.diff = 1e-10;
            param.eIter = 5;
            param.eIter_num = 0.95;
        end
        %data  [A,rows,cols,entries,rep,field,symm]=mmread('illc1033.mtx');
%         \hline 
%  \multirow{4}{*}{$ 1033\times 320 $}& IFM & 7.14143 & 7.72974e-15 & (1500,0)  & 0.12 \\
% & DHA & 7.14143 & 3.53e-16 & (336,1)  & 0.051 \\
% & CHA & 7.14143 & 1.42168e-15 & (12,2)  & 0.018 \\
% & RHA & 7.14143 & 1.42168e-15 & (12,2)  & 0.018 \\
% & PHA & 7.14143 & 2.4265e-15 & (15,2)  & 0.018 \\
         function param = parameter_illc1033()
            param.tol = 1e-15;
            param.maxIter = 600;
            param.nf = 3;
            param.steplengthOrk = 5;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.7;
            param.diff = 1e-3;
            param.eIter = 2;
            param.eIter_num = 0.95;
        end
        %data [A,rows,cols,entries,rep,field,symm]=mmread('illc1850.mtx');
%         \hline 
%  \multirow{4}{*}{$ 1850\times 712 $}& DHA & 9.59166 & 6.02927e-16 & (640,1)  & 0.107 \\
% & CHA & 9.59166 & 6.61387e-14 & (10,1)  & 0.073 \\
% & RHA & 9.59166 & 1.07131e-15 & (50,3)  & 0.096 \\
% & PHA & 9.59166 & 1.05423e-10 & (40,3)  & 0.083 \\
         function param = parameter_illc()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.5;
            param.diff = 1e-5;
            param.eIter = 2;
            param.eIter_num = 0.95;
        end

%          \multirow{4}{*}{$ 2000\times 1400 $}& DHA & 6.57294e-13 & 4.57e-12 & (361,0)  & 5.4334 \\
% & CHA & 2.63079e-14 & 5.18038e-13 & (23,1)  & 1.4294 \\
% & RHA & 6.5965e-13 & 5.91268e-12 & (123,2)  & 3.0509 \\
% & PHA & 2.71702e-14 & 5.36049e-13 & (10,1)  & 1.5809 \\
         function param = parameter2000_1400_10()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.8;
            param.diff = 10*eps;
            param.eIter = 2;
            param.eIter_num = 0.95;
        end
%          \multirow{4}{*}{$ 1000\times 800 $}& DHA & 9.12742e-16 & 1.39023e-14 & (949,0)  & 1.3552 \\
% & CHA & 3.12238e-15 & 4.51886e-14 & (61,3)  & 0.2296 \\
% & RHA & 3.54412e-15 & 5.02664e-14 & (244,4)  & 0.4105 \\
% & PHA & 3.39384e-15 & 4.85255e-14 & (51,4)  & 0.2694 \\
        function param = parameter1000_800_10()
            param.tol = 1e-15;
            param.maxIter = 300;
            param.nf = 3;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.6;
            param.diff = 10*eps;
            param.eIter = 2;
            param.eIter_num = 0.9;
        end
%          \multirow{4}{*}{$ 1000\times 700 $}& DHA & 8.91287e-13 & 4.27203e-12 & (344,0)  & 1.1726 \\
% & CHA & 1.26661e-14 & 1.76672e-13 & (23,1)  & 0.1881 \\
% & RHA & 4.12112e-13 & 2.47838e-12 & (68,2)  & 0.4121 \\
% & PHA & 2.80832e-13 & 1.43003e-12 & (173,1)  & 0.7279 \\
        function param = parameter1000_700_10()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.8;
            param.diff = 10*eps;
            param.eIter = 1;
        end
%          \multirow{4}{*}{$ 1000\times 600 $}& DHA & 1.63235e-14 & 2.09107e-13 & (400,1)  & 1.1151 \\
% & CHA & 1.64732e-14 & 2.14597e-13 & (32,1)  & 0.2979 \\
% & RHA & 1.60181e-14 & 2.06093e-13 & (81,1)  & 0.3486 \\
% & PHA & 1.64198e-14 & 2.09168e-13 & (118,1)  & 0.4843 \\
        function param = parameter1000_600_10()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.9;
            param.diff = 10*eps;
            param.eIter = 2;
        end
        %          \multirow{4}{*}{$ 2000\times 1200 $}& DHA & 2.32273e-15 & 4.50188e-14 & (2400,2)  & 30.507 \\
        % & CHA & 4.51891e-15 & 8.13075e-14 & (310,4)  & 5.671 \\
        % & RHA & 3.91853e-15 & 7.19521e-14 & (320,25)  & 5.327 \\
        % & PHA & 5.16507e-15 & 8.99295e-14 & (370,23)  & 6.343 \\
         function param = parameter200012001()
            param.tol = 1e-15;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.9;
            param.diff = 10*eps;
            param.eIter = 5;
         end
%         \hline 
%  \multirow{4}{*}{$ 2000\times 1200 $}& DHA & 4.79268e-15 & 8.85057e-14 & (1600,2)  & 20.303 \\
% & CHA & 3.86788e-15 & 7.03182e-14 & (900,14)  & 14.557 \\
% & RHA & 4.38813e-15 & 7.85561e-14 & (1000,74)  & 13.927 \\
% & PHA & 4.38813e-15 & 7.85561e-14 & (1000,74)  & 15.246 \\
%    load('test.mat','A','b','x0');
        function param = parameterTest()
            param.tol = 1e-15;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.9;
            param.diff = 10*eps;
            param.eIter = 2;
        end
        % count=10
%          \multirow{4}{*}{$ 2000\times 1200 $}& DHA & 3.49834e-13 & 1.67807e-12 & (935,1)  & 12.5539 \\
% & CHA & 3.12937e-14 & 5.67979e-13 & (30,1)  & 2.1504 \\
% & RHA & 3.05403e-14 & 5.48779e-13 & (105,1)  & 2.2535 \\
% & PHA & 2.14908e-13 & 1.0948e-12 & (403,1)  & 6.2704 \\
         function param = parameter20001200()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.9;
            param.diff = 10*eps;
            param.eIter = 2;
         end
 
  
    end
end