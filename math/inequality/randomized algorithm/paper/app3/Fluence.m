function [x]=Fluence(A,B,l,u,k,r)
[m1,n]=size(A);
[m2,n]=size(B);
H=[1/sqrt(m2)*B;sqrt(r/m1)*A];
w=0;
dw=[1/sqrt(m2)*l;sqrt(r/m1)*(u+w)];
              f = -H'*dw;
              HTH=H'*H;
            if strcmp(prob.nnls,'quadprog')   
                options = optimoptions(@quadprog,'Display','off');
                x = quadprog(HTH,f,[],[],[],[],lb,ub,[],options);
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
            
           
%             step = prob.structs{ii}.terms{jj}.step;
%             nVoxels = prob.structs{ii}.nVoxels;
%             coeff = step*weight/nVoxels;
            coeff = 1;
            % Calculate gradient step
%             dose = prob.structs{ii}.A*prob.x;
            %res = (-1)^s*(H*x - d);
            res=H*x - d;
            wPrev = w;
            wStep = wPrev + coeff*(res - wPrev);
            
            % Project onto set ||(w)_+||_0 <= k
            wProj = projW(wStep,k);
%             wDiff = norm(wProj - wPrev)/step;
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