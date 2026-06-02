function [center, radius] = circleFrom3Points(p1, p2, p3)
    % Compute the circle defined by three 2-D points.
    mid1 = (p1 + p2) / 2;
    mid2 = (p2 + p3) / 2;

    dir1 = p2 - p1;
    dir2 = p3 - p2;

    perp1 = [-dir1(2), dir1(1)];
    perp2 = [-dir2(2), dir2(1)];

    A = [perp1' -perp2'];
    b = (mid2 - mid1)';

    ts = A\b;
    t = ts(1);

    center = mid1 + t * perp1;
    radius = norm(center - p1);
end
