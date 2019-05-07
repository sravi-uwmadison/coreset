function [gu, gv] = numerical_gradient(A, B, u, v, h)
    if ~exist('h', 'var')
        h = 1e-6;
    end

    obj = @(u, v) norm(A * u - B * v, Inf);
    
    
    m = length(u);
    gu = zeros(m, 1);
    d = zeros(m,1);
    for i=1:m
        d(i) = 1;
        gu(i) = obj(u + h*d, v) - obj(u - h*d, v);
        d(i) = 0;
    end
    gu = gu / (2*h);
    
    n = length(v);
    gv = zeros(n, 1);
    d = zeros(n,1);
    for i=1:n
        d(i) = 1;
        gv(i) = obj(u, v + h * d) - obj(u, v - h * d);
        d(i) = 0;
    end
    gv = gv / (2*h);
end