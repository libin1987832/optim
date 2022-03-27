        
f = -A'*d;
            if strcmp(prob.nnls,'quadprog')   
                options = optimoptions(@quadprog,'Display','off');
                x = quadprog(H,f,[],[],[],[],lb,ub,[],options);
            else
                func = @(x)FluenceMapOpt.quadObj(x,H,f);
                options.verbose = 0;
                options.method = 'newton';
                x = minConf_TMP(func,x0,lb,ub,options);
            end
        end

        
        
        function wDiff = updateW(prob,ii,jj)
            % UPDATEW Update proxy variable.
            
            % Grab variables
            s = strcmp(prob.structs{ii}.terms{jj}.type,'ldvc');
            k = prob.structs{ii}.terms{jj}.k;
            step = prob.structs{ii}.terms{jj}.step;
            weight = prob.structs{ii}.terms{jj}.weight;
            nVoxels = prob.structs{ii}.nVoxels;
            coeff = step*weight/nVoxels;
            
            % Calculate gradient step
            dose = prob.structs{ii}.A*prob.x;
            res = (-1)^s*(dose - prob.structs{ii}.terms{jj}.d);
            wPrev = prob.structs{ii}.terms{jj}.w;
            wStep = wPrev + coeff*(res - wPrev);
            
            % Project onto set ||(w)_+||_0 <= k
            wProj = FluenceMapOpt.projW(wStep,k);
            wDiff = norm(wProj - wPrev)/step;
            prob.structs{ii}.terms{jj}.w = wProj;
        end