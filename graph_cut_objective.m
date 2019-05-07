function S = graph_cut_objective(x, W)
    % S = graph_cut_objective(x, W)

    [Wi, Wj, Wv] = find(W);
    E = length(Wi);
    
    S = 0;
    for i=1:E
        S = S + Wv(i) * norm(x(:, Wi(i)) - x(:, Wj(i)), 1);
    end
end