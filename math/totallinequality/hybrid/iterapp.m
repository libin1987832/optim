function y = iterapp(op,afun,atype,afcnstr,x,varargin)
if strcmp(atype,'matrix')
    switch lower(op)
        case 'mtimes'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' * x;
            else
                y = afun * x;
            end
        case 'mldivide'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' \ x;
            else
                y = afun \ x;
            end
        otherwise
            error(message('MATLAB:iterapp:InvalidOp'))
    end
else
    try
        if (nargin >= 6) && isequal(varargin{end},'notransp')
            % A request for A*x coming from BICG, LSQR and QMR
            try
                % New syntax: we now request afun(x,P1,P2,...,PN,'notransp')
                y = afun(x,varargin{:});
            catch
                % Old syntax: we used to accept afun(x,P1,P2,...,PN)
                y = afun(x,varargin{1:end-1});
            end
        else
            % A request for A*x
            % coming from BICGSTAB, CGS, GMRES, MINRES, PCG or SYMMLQ
            % with the call always afun(P1,P2,...,PN)
            % or a request for A'*x coming from
            % BICG, LSQR and QMR in the afun(x,P1,P2,...,PN,'transp') case
            y = afun(x,varargin{:});
        end
    catch ME
        error(message('MATLAB:iterapp:InvalidInput', atype,afcnstr, ME.message));
    end

    if ~iscolumn(y)
        error(message('MATLAB:iterapp:MustReturnColumn', atype, afcnstr));
    end
end
