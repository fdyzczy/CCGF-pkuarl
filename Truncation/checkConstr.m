function isViol = checkConstr(p, constr)
% Return a 1-by-N indicator for samples that fall in the inadmissible set.

    N = size(p, 2);
    isViol = false(1, N);

    if isempty(constr)
        isViol = double(isViol);
        return;
    end

    if ~iscell(constr)
        constr = {constr};
    end

    for i = 1:length(constr)
        if constr{i}.type == "cur"
            c = constr{i}.center(:);
            r = constr{i}.a;
            newp = p(1:2, :);
            dist = sqrt(sum((newp - c).^2, 1));
            isViol = isViol | (dist <= r);

        elseif constr{i}.type == "lin"
            a = constr{i}.a(:);
            b = constr{i}.b;
            isViol = isViol | (a' * p <= b);

        elseif constr{i}.type == "pol"
            newp = p(1:2, :);
            vtList = constr{i}.vtList;
            inPol = inpolygon(newp(1, :), newp(2, :), vtList(:, 1), vtList(:, 2));
            isViol = isViol | inPol;
        end
    end

    isViol = double(isViol);
end
