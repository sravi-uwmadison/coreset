function D = diaminf(A)
    % D = diaminf(A)
    %
    %   Finds the infinity-norm diameter of the columns of A
    %
    
    D = max(max(A, [], 2) - min(A, [], 2));
end