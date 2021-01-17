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
            param.eIter = 10;
        end
    end
end