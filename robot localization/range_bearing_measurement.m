function [z, bounds] = range_bearing_measurement(x, lm)
    % x is 2-by-N; each column is a robot position.
    % lm is a single landmark as 1-by-2 or 2-by-1.
    if size(lm,1) == 1
        lm = lm';
    end
    dx = lm(1) - x(1,:);
    dy = lm(2) - x(2,:);
    range = sqrt(dx.^2 + dy.^2);
    bearing = atan2(dy, dx);
    z = [range; wrapToPi(bearing)];

    bounds = [ -Inf, Inf;
               -pi,  pi ];
end
