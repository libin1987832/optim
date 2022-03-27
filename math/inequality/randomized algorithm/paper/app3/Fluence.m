        
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

        
        
        function wDiff = updateW(w,k,step,weight)
            % UPDATEW Update proxy variable.
            
            % Grab variables
            
           
            step = prob.structs{ii}.terms{jj}.step;
            nVoxels = prob.structs{ii}.nVoxels;
            coeff = step*weight/nVoxels;
            
            % Calculate gradient step
%             dose = prob.structs{ii}.A*prob.x;
            %res = (-1)^s*(H*x - d);
            res=H*x - d;
            wPrev = w;
            wStep = wPrev + coeff*(res - wPrev);
            
            % Project onto set ||(w)_+||_0 <= k
            wProj = projW(wStep,k);
            wDiff = norm(wProj - wPrev)/step;
            w = wProj;
        end
         function w = projW(w,k)
            % PROJW Project w onto the set satisfying ||max(0,w)||_0 <= k.
            idxPos = w > 0;
            if sum(idxPos) > k
                wPos = w(idxPos);
                [~,idxSort] = sort(wPos,'descend');
                wPos(idxSort(k+1:end)) = 0;
                w(idxPos) = wPos;
            end
        end