function alpha = searchstep(tao, z, u)
    alpha  = 0;
    if z >= 1
        alpha = 1;
    else
        alpha = z;
    end
    if tao > 0 && z < 1 && z >u
        alpha = u;
    end
%    if tao>0
%             if z>=1
%                alpha=1;
%             elseif z<1&&z>u
%                     alpha=z;
%             else
%                 alpha=u;
%             end
%             
%    else
%             if z>=1
%                 alpha=1;
%             else
%                 alpha=z;
%             end
%     end