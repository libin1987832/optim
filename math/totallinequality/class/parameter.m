classdef parameter
    % fixedMprgp
    properties
        fixed_tol
        fixed_maxit
    end
    % MPRGP
    properties
        mprgp_L
        mprgp_a
        mprgp_delta
        mprgp_Ftol
        mprgp_maxIter
    end
    methods
         function param = parameter()
            param.fixed_tol = 1e-10;
            param.fixed_maxit = 1000;
            param.mprgp_L = 1;
            param.mprgp_a = 0.01;
            param.mprgp_delta = 1e-30;
            param.mprgp_Ftol = 1e-15;
            param.mprgp_maxIter = 100;
         end
    end
end