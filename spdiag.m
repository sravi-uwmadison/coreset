function A = spdiag(d)
    % A = spdiag(d)

    n = length(d);
    A = sparse(1:n, 1:n, d);
end