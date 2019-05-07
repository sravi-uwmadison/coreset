function median_subproblem(x, dfS, szT, n)
    obj = zeros(n+1, n+1);
    
    vmin = Inf;
    vmax = -Inf;
    
    zmin = x;
    
    for i=1:(n+1)
        for j=1:(n+1)
            z = zeros(3, 1);
            
            z(1) = (i-1)/n;
            z(2) = (j-1)/n;
            if z(1) + z(2) > 1
                obj(i,j) = -Inf;
                continue
            end
            z(3) = 1 - z(1) - z(2);
            
            
            obj(i,j) = dfS' * z + szT * norm(z - x);
            
            if obj(i, j) < vmin
                vmin = obj(i,j);
                zmin = z;
            end
            vmax = max(vmax, obj(i,j));
        end
    end
    
    zmin
    contourf((0:n)/n, (0:n)/n, obj, vmin:((vmax - vmin)/50):vmax);
    
    % Produces zmin in the interior for
    % median_subproblem([1;0;1], [0; -2; -0]/sqrt(5), 1, 100)
    
    % Another ugly case:
    % median_subproblem([1;0;0], [0; -1; -1]*3, 2, 100)
end