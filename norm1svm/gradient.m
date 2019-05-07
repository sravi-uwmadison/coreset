function [gu, gv] = gradient(A, B, u, v)
    diff = A*u - B*v;
    
    nrm = norm(diff, Inf);
    
    dgrad = zeros(length(diff), 1);
    dgrad(diff == nrm) = 1;
    dgrad(diff == -nrm) = -1;
    
    gu = A' * dgrad;
    gv = -B' * dgrad;
end