function [E, L, Df] = converge_bounds(X, y)
    Df = 2 * (diaminf(X(y < 0)) + diaminf(X(y >= 0)))^2;
    L = norm(full(X), 2);
    E = 2^(5/2) * L + 2^(3/2) * Df;
end