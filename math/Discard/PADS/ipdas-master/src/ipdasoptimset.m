function option = ipdasoptimset(varargin)
%===========================================================================
% Description	: Option structure for ipdas
% Parameters	: method	| ssm solver options
%                           |- 'exact'
%                           |- 'bound1'
%                           |- 'bound2'(default)
%                           |- 'inexact'
%               : tol_opt 	| tolerance for kkt residual(default 1e-10)
%               : tol_iter	| maximum number of iterations(default 1e+3)
%               : tol_res_i	| tolerance for kkt residual on I set(default 1e-10)
%               : display   | display iterations or not
%                           |- 'on'(default)
%                           |- 'off'
% Author        : Zheng Han
% Date          : 12/12/2012
% Example       : option = ipdasoptimset('method','inexact','tol_iter',1000)
%===========================================================================

p = inputParser;
p.addParamValue('method','bound2',@(x)any(strcmpi(x,{'exact','bound1','bound2','inexact'})));
p.addParamValue('tol_opt',1e-10,@isnumeric);
p.addParamValue('tol_iter',1e+2,@isnumeric);
p.addParamValue('tol_res_i',1e-10,@isnumeric);
p.addParamValue('display','on',@(x)any(strcmpi(x,{'on','off'})));
p.parse(varargin{:});
option = p.Results;

end