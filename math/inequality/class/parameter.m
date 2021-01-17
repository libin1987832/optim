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
    end
    methods
        function param = parameter()
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
        %          \multirow{4}{*}{$ 2000\times 1200 $}& DHA & 2.32273e-15 & 4.50188e-14 & (2400,2)  & 30.507 \\
        % & CHA & 4.51891e-15 & 8.13075e-14 & (310,4)  & 5.671 \\
        % & RHA & 3.91853e-15 & 7.19521e-14 & (320,25)  & 5.327 \\
        % & PHA & 5.16507e-15 & 8.99295e-14 & (370,23)  & 6.343 \\
         function param = parameter20001200()
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
    end
end