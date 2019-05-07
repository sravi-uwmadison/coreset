function [w, b] = norm1svm_primal_from_dual(X, y, alpha)
    % [w, b] = norm1svm_primal_from_dual(X, y, alpha)

    % TODO: non-seperable case
    
    [d, n] = size(X);
    
    sv = abs(alpha) > 1e-8;
    wb = [X(:, sv)' ones(sum(sv),1)] \ y(sv);
    w = wb(1:d);
    b = wb(end);
end