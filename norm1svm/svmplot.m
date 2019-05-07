function svmplot(X, y, w, b)
    assert(size(X, 1) >= 2);
    assert(size(X, 2) == length(y));
    assert(size(X, 1) == length(w));

    hold on
    scatter(X(1, y < 0), X(2, y < 0), 'r+');
    scatter(X(1, y >= 0), X(2, y >= 0), 'bo');
    
    square_axes
    
    % Plot decision boundary
    v = axis;
    xmin = v(1); xmax = v(2); ymin = v(3); ymax = v(4);
    
    if abs(w(1)) > abs(w(2))
        x1 = -(w(2) * ymin + b)/w(1);
        x2 = -(w(2) * ymax + b)/w(1);
        plot([x1, x2], [ymin, ymax], 'k-');
       
    else
        y1 = -(w(1) * xmin + b)/w(2);
        y2 = -(w(1) * xmax + b)/w(2);
        plot([xmin, xmax], [y1, y2], 'k-');
    end
    hold off
end