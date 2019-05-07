function d = svmT(x, eps)
    n = length(x);
    
    d = [];
    
    for i=1:n
        u = x - sign(x) * eps;
        
        u(i) = x(i) - eps;
        d = [d; infnormsubdiff(u)];
        
        u(i) = x(i) + eps;
        d = [d; infnormsubdiff(u)];
    end
    
    d = unique(d, 'rows');
end