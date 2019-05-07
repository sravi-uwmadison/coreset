function a = binary_line_search(x, d, f, amin, amax, eps)
    if ~exist('eps', 'var')
        eps = 1e-6;
    end
    
    epsilon = (amax - amin) * eps;
    while amax - amin > epsilon
        amid = (amax + amin)/2;
        
        if f(x + d*amax) < f(x + d*amin)
            amin = amid;
        else
            amax = amid;
        end
    end
    
    a = (amax + amin)/2;
end