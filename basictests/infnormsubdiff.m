function d = infnormsubdiff(u)
    d = zeros(size(u));
    A = abs(u) == norm(u, Inf);
    
    d(A) = sign(u(A));
end