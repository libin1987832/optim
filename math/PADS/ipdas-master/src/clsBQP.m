classdef clsBQP < handle
    %==================================================
    % Description	: Class for bound constrained QP of the form
    %               : min 0.5x'Hx + c'x
    %               : st     x <= ub
    % Algorithm     : Designed for primal-dual-active-set
    % Author        : Zheng Han
    % Date          : 12/09/2012
    %==================================================
    properties(SetAccess = private, Hidden = true)
        H;
        c;
        ub;
    end
    %--------------------------------------------------
    properties(SetAccess = public)
        x;           % primal solution
        y;           % dual solution
        U;           % upper-active-set
        I;           % inactive-set
        cg;          % number of CG sub-iterations
        iter;        % number of pdas-iterations
        status;      % exit state: optimal(1), not-optimal(0), exceed-limit(-1)
    end
    %--------------------------------------------------
    properties(Dependent = true, SetAccess = private)
        n;           % problem size
        f;           % objective value
        kkt;         % kkt residual
        invH;        % inverse of H
        res_i;       % kkt residual of I set part
    end
    %==================================================de
    methods
        %--------------------------------------------------
        function obj = clsBQP(H,c,UB,U,I)
            obj.H = H;
            obj.c = c;
            obj.ub = UB;
            obj.U = U;
            obj.I = I;
            obj.cg = 0;
            obj.iter = 0;    
            obj.status = 0;
            if isempty(obj.x)
                obj.x = zeros(length(c),1);
            end
            if isempty(obj.y)
                obj.y = -(obj.H*obj.x + obj.c);
            end
            if isempty(obj.I)
                obj.I = zeros(1,0);                
            end
            if isempty(obj.U)
                obj.U = zeros(1,0);
            end
        end
        %--------------------------------------------------
	% Class copy constructor, refer
	% http://www.mathworks.com/matlabcentral/newsreader/view_thread/257925
	function new = copy(this)
	    % Instantiate new object of the same class
	    new = feval(class(this),this.H,this.c,this.ub,this.U,this.I);

	    % Copy all non-hidden properties
%	    p = properties(this);
%	    for i = 1:length(p)
%	        new.(p{i}) = this.(p{i});
%	    end
	end
        %--------------------------------------------------
        function n = get.n(obj)
            n = length(obj.c);
        end
        function f = get.f(obj)
            f = 0.5*obj.x'*obj.H*obj.x + obj.c'*obj.x;
        end
        function invH = get.invH(obj)
            invH = inv(obj.H);
        end
        function res_i = get.res_i(obj)
            if ~isempty(obj.I)
                res_i = obj.H(obj.I,:)*obj.x + obj.c(obj.I);
            else
                res_i = 0;
            end
        end
        function kkt = get.kkt(obj)
            kkt = max(norm(obj.res_i,inf),norm( min([obj.y;obj.ub - obj.x]) ,inf));
        end

        %--------------------------------------------------
        function obj = updateVars(obj,method,tol_res_i)
            % Set fixed variables
            obj.x(obj.U) = obj.ub(obj.U);
            obj.y(obj.I) = 0;
            
            % Compute x(I) and y(U) (exactly/inexactly)
            if strcmpi(method,'exact')
                % Compute exact solution
                if ~isempty(obj.I)
                    obj.x(obj.I) = -obj.H(obj.I,obj.I)\(obj.H(obj.I,obj.U)*obj.x(obj.U) + obj.c(obj.I));
                end
                if ~isempty(obj.U)
                    obj.y(obj.U) = -(obj.H(obj.U,:)*obj.x + obj.c(obj.U));
                end

            else
                % Inexact methods with guarantees
                obj.cg = 0;
                r0 = obj.res_i;
                p = -obj.res_i;
                while norm(obj.res_i,inf) > tol_res_i
                    % Do a CG iteration
                    a = obj.res_i'*obj.res_i/(p'*obj.H(obj.I,obj.I)*p);
                    obj.x(obj.I) = obj.x(obj.I) + a*p;
                    obj.y(obj.U) = -(obj.H(obj.U,:)*obj.x + obj.c(obj.U));
                    r1 = r0 + a*obj.H(obj.I,obj.I)*p;
                    obj.cg = obj.cg + 1;
                    % Check if the solution is accurate enough
                    if strcmpi(method,'bound1')
                        if (isempty(obj.I) || norm(obj.res_i,inf) <= min(abs(obj.ub(obj.I)-obj.x(obj.I)))/norm(obj.invH(obj.I,obj.I),inf) )...
                                && (isempty(obj.U) || (norm(obj.res_i,inf) <= min(abs(obj.y(obj.U)) )/(norm(obj.H(obj.U,obj.I),inf)*norm(obj.invH(obj.I,obj.I),inf)) ))
                            break;
                        end
                    elseif strcmpi(method,'bound2')
                        s = norm(obj.H(obj.I,obj.I),inf)+1;
                        coef = s*(1-norm(eye(length(obj.I)) - obj.H(obj.I,obj.I)/s,inf));
                        
                        if (isempty(obj.I) || norm(obj.res_i,inf) <= coef*min( abs(obj.ub(obj.I) - obj.x(obj.I)) ) )...
                                && (isempty(obj.U) || norm(obj.res_i,inf) <= coef*min( abs(obj.y(obj.U)) )/norm(obj.H(obj.U,obj.I),inf))
                            break;
                        end
                    end
                    
                    % Prepare for next CG iteration
                    b = r1'*r1/(r0'*r0);
                    r0 = r1;
                    p = -r1 + b*p;

                end
                
            end
        end
        %--------------------------------------------------
        function obj = updatePartition(obj)
            % Violated indices
            VP = obj.I(  obj.x(obj.I) > obj.ub(obj.I)  );
            VD = obj.U(  obj.y(obj.U) < 0 );
            
            % Reset the index set
            obj.I = union(VD,setdiff(obj.I,VP));
            obj.U = union(VP,setdiff(obj.U,VD));
            obj.iter = obj.iter + 1;
        end
        %--------------------------------------------------
        function print(obj)
            if obj.iter == 0
                fprintf('============================================================\n');
                fprintf('    k      obj         kkt        |I|    |U|   |CG|   cg_res\n');
                fprintf('============================================================\n');
            end
            fprintf('%5d  %+.4e  %.4e  %5d  %5d  %5d  %5d  %.4e\n',...
                obj.iter,obj.f,obj.kkt,length(obj.I),length(obj.U),obj.cg,norm(obj.res_i,inf));
            fprintf('\n');
        end
    end
    %==================================================
end