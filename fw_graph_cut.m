function [x obj xs] = fw_graph_cut(u, W, T)
tic
    sizeW = size(W,1);
    sizeT = length(T);
    % x = fw_graph_cut(u, W, x0)
    %
    %   Solve a graph-cut problem with unary weights u and weighted
    %   adjacency matrix W and terminals T.
    %
    
    m = 100; % Restricting the number of terminals to keep it simple        
    start = 50;
    % Test case
%     try
%         T = T(9000:9168,:); % Choose the first m terminals    
%     end
%     try
%         T = T(:,9000:9168); % Choose the m terminals    
%     end
    n = size(W, 1);    
    d = length(T);
%     if ~exist(x0, 'var')
%         x0 = zeros(n, 1);
%     end
    x0 = 0*ones(n*d, 1); % Start from the origin though this is not needed
    % Simplex constraints for all the nodes
    % Equality matrix     
     i = reshape(repmat([1:n],d,1),n*d,1);
     j = [1:n*d]';
     Aeq = sparse(i,j,1);
     
    % Decompose equality matrix into two inequality matrices and add the
    % nonnegativity constraints
     
    % Create the sparse identity of size n*d and d
     s_eye_n_d = sparse([1:n*d],[1:n*d],1);
     s_eye_d = sparse([1:d],[1:d],1);    
     A = [Aeq;-Aeq;-s_eye_n_d];
     b = [ones(size(Aeq,1),1);-ones(size(Aeq,1),1) ; sparse(n*d,1)];
     
     % Impose that the terminals T are basis vectors of d dimensional
     % vector space or in other words the vertices of the d dimensional
     % simplex for a general T. The ordering of nodes is preserved.
    ind_row = [1:d*d];    
    ind_col = zeros(size(ind_row,1),1);
    for k=1:d
        ind_col((k-1)*d+1:k*d)=[(T(k)-1)*d+1:T(k)*d];
    end
    basis_cons_eq = sparse(ind_row,ind_col,1);    
    % Add zeros if necessary to make the matrix consistent with the size
    if size(A,2)-size(basis_cons_eq,2) > 0    
        zero_mat = sparse(d*d,size(A,2)-size(basis_cons_eq,2));    
        basis_cons_eq = [basis_cons_eq, zero_mat];
    end
    basis_cons = sparse([basis_cons_eq; -basis_cons_eq]);    
    basis = reshape(s_eye_d,d*d,1);    
    A = [A;basis_cons];
    b = [b;basis;-basis];      
    clear basis basis_cons i j basis_cons_eq zero_mat ind_col s_eye_n_d
    clear s_eye_d ind_row
    % Basis constraints for terminals
    u = zeros(size(A,2),1);    
    W = W - diag(diag(W));
    df = @(x) [u, u] + graph_cut_subdifferential(W, x);      
    [x] = boxsubgradient(A, b, df, x0, 0, 10);
    toc
end