function alpha = searchstep(tao, z, u)
    if z >= 1
        alpha = 1;
    else
        alpha = z;
    end
    if tao > 0 && z < 1 && z <u  
        alpha = u;
    end
    
end
