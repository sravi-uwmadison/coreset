function M = incidence_matrix(W)
    % M = incidence_matrix(W)
    %
    %   Return a weighted incidence matrix from weight adjacency matrix W
    %
    
    n = size(W, 1);
    [Wi, Wj, Wv] = find(W);
    E = length(Wi);
    
    Mi = reshape([1:E; 1:E], 2*E, 1);
    Mj = reshape([Wi'; Wj'], 2*E, 1);
    Mv = reshape([Wv'; -Wv'], 2*E, 1);
    
    M = sparse(Mi, Mj, Mv, E, n);
end