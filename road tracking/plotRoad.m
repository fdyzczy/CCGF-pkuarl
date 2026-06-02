function roadhdl = plotRoad(avgPoints)

    hold on
    arc_indices = {
        [12,14,2], [12,25,32], [32,42,31], [31,20,2], ...
        [9,7,6], [6,22,30], [30,40,36], [36,17,9]
    };
    for g = 1:length(arc_indices)
        pts = avgPoints(arc_indices{g},:);
        roadhdl = plotArcThrough3Points(pts(1,:), pts(2,:), pts(3,:), 'g', 'LineWidth', 2);
    end

    axis equal;
end


function roadhdl =plotArcThrough3Points(ptA, ptB, ptC, varargin)
    pts = [ptA; ptB; ptC];

   [center, r] = circleFrom3Points(ptA, ptB, ptC);

    roadhdl =plotShortArc(center, r, pts(1,:), pts(2,:), 'k', 'LineWidth', 2);
    roadhdl =plotShortArc(center, r, pts(2,:), pts(3,:), 'k', 'LineWidth', 2);
end


function roadhdl =plotShortArc(center, r, pa, pb, varargin)
    theta_a = atan2(pa(2) - center(2), pa(1) - center(1));
    theta_b = atan2(pb(2) - center(2), pb(1) - center(1));
    if theta_b<theta_a
        theta_b = theta_b+2*pi;
    end
    dtheta = mod(theta_b - theta_a, 2*pi);
    if dtheta > pi
        theta_range = linspace(theta_a, theta_b - 2*pi, 100);
    else


        theta_range = linspace(theta_a, theta_b, 100);
    end
    xx = center(1) + r * cos(theta_range);
    yy = center(2) + r * sin(theta_range);
    roadhdl = plot(xx, yy, varargin{:},"Displayname","");
end
