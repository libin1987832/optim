classdef parametern
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
         function param = parametern()
            param.fixed_tol = 1e-13;
            param.fixed_maxit = 2000;
            param.mprgp_L = 1;
            param.mprgp_a = 0.5;
            param.mprgp_delta = 1e-30;
            param.mprgp_Ftol = 1e-25;
            param.mprgp_maxIter =20;
         end
         %param.mprgp_a = 1/norm(A'*A, inf)*10;
%         [A,b,x0] = readData(7,500,500,300);
%         & FMprgp & 18.19 & 0 & 9.15974e-14 & 8.31788e-13 & 0.498  \\
% & ls & 18.19 & 5.49046e-22 & 1.11788e-13 & -1.36082e-12 &0.643137 \\
         function param = parametern1()
            param.fixed_tol = 1e-13;
            param.fixed_maxit = 1000;
            param.mprgp_L = 1;
            param.mprgp_a = 0.5;
            param.mprgp_delta = 1e-30;
            param.mprgp_Ftol = 1e-25;
            param.mprgp_maxIter = 8;
         end
        %         & FMprgp & 60.2694 & 0 & 9.42747e-12 & -2.64471e-12 & 0.49  
% & ls & 60.2694 & 2.6693e-12 & 3.94438e-05 & -7.83655e-07 &0.254098 \\
% & IPG & 60.2694 & 0 & 2.97908e-11 & -5.64708e-12 &0.976 \\
% save('test_paper2.mat','A','b','x0')
         function param = parameter2()
            param.fixed_tol = 1e-11;
            param.fixed_maxit = 1000;
            param.mprgp_L = 1;
            param.mprgp_a = 1;
            param.mprgp_delta = 1e-50;
            param.mprgp_Ftol = 1e-10;
            param.mprgp_maxIter = 6;
         end

        %         & FMprgp & 66.1538 & 0 & 9.94362e-12 & -1.06281e-12 & 1.385  
% & ls & 66.1538 & 5.23672e-14 & 1.11883e-06 & -1.3101e-12 &0.247618 \\
% & IPG & 66.1538 & 0 & 3.78163e-11 & -2.20805e-11 &5.002 \\
% save('test_paper1.mat','A','b','x0')
         function param = parameter1()
            param.fixed_tol = 1e-11;
            param.fixed_maxit = 1000;
            param.mprgp_L = 1;
            param.mprgp_a = 0.5;
            param.mprgp_delta = 1e-50;
            param.mprgp_Ftol = 1e-10;
            param.mprgp_maxIter = 10;
         end
    end
end