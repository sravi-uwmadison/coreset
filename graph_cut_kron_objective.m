function S = graph_cut_kron_objective(x, W)
    % graph_cut_kron_objective(x, W)

    [d n] = size(x);
    B = incidence_matrix(W);
    E = size(B, 1);
    %S = norm(kron(B, eye(d)) * reshape(x, d*n, 1), 1);
    S = norm(reshape(x*B', n * E, 1), 1);
end