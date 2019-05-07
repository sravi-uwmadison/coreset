function Sduff = graph_cut_objective_diff(x1, x2, W)
    % Sdiff = graph_cut_objective_diff(x1, x2, W)
    %
    %   Returns f(x1) - f(x2), where f is the graph cut objective  
    %

    [Wi, Wj, Wv] = find(W);
    E = length(Wi);
    
    Sdiff = 0
    for i=1:E
        Si1 = Wv(i) * norm(x1(:, Wi(i)) - x1(:, Wj(i)), 1);
	Si2 = Wv(i) * norm(x2(:, Wi(i)) - x2(:, Wj(i)), 1);
        Sdiff = Sdiff + (Si1 - Si2);
    end
end