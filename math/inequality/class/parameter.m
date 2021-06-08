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
%         \hline 
%  \multirow{4}{*}{$ 2000\times 1600 $}& DHA & 5.15314e-13 & 5.02499e-12 & (186,0)  & 4.2586 \\
% & DHA($ \mu=n_f $) & 1.55017e-14 & 3.29507e-13 & (16,2)  & 2.853 \\
% & CHA & 1.55017e-14 & 3.29507e-13 & (16,2)  & 2.852 \\
% & RHA & 1.40002e-14 & 2.9747e-13 & (49,1)  & 1.583 \\
% & PHA & 1.53412e-14 & 3.24847e-13 & (17,2)  & 2.6356 \\
         function param = parameter2000_1600()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 3;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 10;
            param.con1 = 0.95;
            param.con2 =0.6;
            param.diff = 1e-10;
            param.eIter = 2;
            param.eIter_num = 0.95;
        end

%          \multirow{4}{*}{$ 2000\times 1400 $}& DHA & 6.57294e-13 & 4.57e-12 & (361,0)  & 5.4334 \\
% & CHA & 2.63079e-14 & 5.18038e-13 & (23,1)  & 1.4294 \\
% & RHA & 6.5965e-13 & 5.91268e-12 & (123,2)  & 3.0509 \\
% & PHA & 2.71702e-14 & 5.36049e-13 & (10,1)  & 1.5809 \\
         function param = parameter2000_1400()
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
% \hline 
%  \multirow{5}{*}{$ 1000\times 800 $}& DHA & 3.01847e-14 & 2.81707e-13 & (146,0)  & 0.3372 \\
% & DHA($ \mu=n_f $) & 4.68391e-14 & 4.00962e-13 & (119,12)  & 0.4732 \\
% & CHA & 3.26384e-14 & 2.98388e-13 & (142,2)  & 0.3524 \\
% & RHA & 3.01573e-14 & 2.76915e-13 & (128,10)  & 0.4462 \\
% & PHA & 3.01847e-14 & 2.81707e-13 & (146,2)  & 0.4052 \\
        function param = parameter1000_800()
            param.tol = 1e-13;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 10;
            param.lsqrIter = 1;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 1;
            param.con1 = 0.95;
            param.con2 =0.7;
            param.diff = eps;
            param.eIter = 10;
            param.eIter_num = 1;
        end
%         \hline 

% \hline 
% \multirow{5}{*}{$ 1000\times 700 $}& DHA & 6.91526e-14 & 4.20527e-13 & (321,0)  & 2.3784 \\
% & DHA($ \mu=n_f $) & 4.29225e-14 & 5.97179e-13 & (107,11)  & 1.7088 \\
% & CHA & 2.30014e-14 & 2.68618e-13 & (120,2)  & 1.2836 \\
% & RHA & 1.48962e-14 & 2.01582e-13 & (112,6)  & 1.2998 \\
% & PHA & 9.4438e-15 & 1.29351e-13 & (148,3)  & 1.4854 \\
        function param = parameter1000_700_3()
            param.tol = 1e-13;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 40;
            param.lsqrIter = 1;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 1;
            param.con1 = 0.95;
            param.con2 =0.838;
            param.diff = eps;
            param.eIter = 5;
            param.eIter_num = 1;
        end
%         \hline 
% \multirow{5}{*}{$ 1000\times 800 $}& DHA & 6.1367e-14 & 5.29394e-13 & (145,0)  & 2.5368 \\
% & DHA($ \mu=n_f $) & 3.91247e-14 & 5.85966e-13 & (49,10)  & 1.4232 \\
% & CHA & 2.17999e-14 & 2.71852e-13 & (81,4)  & 1.2486 \\
% & RHA & 2.25086e-14 & 3.46644e-13 & (63,6)  & 1.0244 \\
% & PHA & 1.65375e-14 & 2.55399e-13 & (71,2)  & 1.1024 \\
                  function param = parameter_1000_800()
            param.tol = 1e-13;
            param.maxIter = 300;
            param.nf = 5;
            param.steplengthOrk = 50;
            param.lsqrIter = 1;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 1;
            param.con1 = 0.95;
            param.con2 =0.6;
            param.diff = eps;
            param.eIter = 2;
            param.eIter_num = 1;
                  end
%                   \hline 
%  \multirow{5}{*}{$ 1000\times 600 $}& DHA & 9.00995e-15 & 1.15379e-13 & (400,1)  & 3.1542 \\
% & DHA($ \mu=n_f $) & 2.93706e-14 & 3.69878e-13 & (38,8)  & 2.6816 \\
% & CHA & 9.60513e-15 & 1.20571e-13 & (83,4)  & 1.844 \\
% & RHA & 1.87502e-14 & 2.43415e-13 & (102,3)  & 1.5406 \\
% & PHA & 1.10833e-14 & 1.39912e-13 & (123,1)  & 1.4712 \\ 30 
           function param = parameter_600()
            param.tol = 1e-13;
            param.maxIter = 300;
            param.nf = 5;
            param.steplengthOrk = 3;
            param.lsqrIter = 1;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 2;
            param.con1 = 0.95;
            param.con2 =0.6;
            param.diff = eps;
            param.eIter = 3;
            param.eIter_num = 1;
          end
           function param = parameter_2()
            param.tol = 1e-13;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 30;
            param.lsqrIter = 1;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 1;
            param.con1 = 0.95;
            param.con2 =0.86;
            param.diff = eps;
            param.eIter = 10;
            param.eIter_num = 1;
        end
%          \multirow{4}{*}{$ 1000\times 700 $}& DHA & 8.91287e-13 & 4.27203e-12 & (344,0)  & 1.1726 \\
% & CHA & 1.26661e-14 & 1.76672e-13 & (23,1)  & 0.1881 \\
% & RHA & 4.12112e-13 & 2.47838e-12 & (68,2)  & 0.4121 \\
% & PHA & 2.80832e-13 & 1.43003e-12 & (173,1)  & 0.7279 \\
        function param = parameter1000_700()
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
        function param = parameter1000_600()
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
         function param = parameterEE()
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
         function param = parameter()
            param.tol = 1e-12;
            param.maxIter = 300;
            param.nf = 10;
            param.steplengthOrk = 20;
            param.lsqrIter = 5;
            param.lsqrTol = 1e-10;
            param.smIsqmIter = 3;
            param.con1 = 0.95;
            param.con2 =0.9;
            param.diff = 10*eps;
            param.eIter = 5;
         end
 
  
    end
end