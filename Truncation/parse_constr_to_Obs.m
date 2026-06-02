function Obs = parse_constr_to_Obs(constr, xdim)
% parse_constr_to_Obs  Convert constr cell-array (pol/lin/cur) to Obs format.
% Output:
%   Obs(i).A : xdim x m
%   Obs(i).b : m x 1
% Interpreted as: Obs(i).A' * x <= Obs(i).b

    if nargin < 2 || isempty(xdim)
        xdim = 2;
    end
    if ~isscalar(xdim) || xdim ~= round(xdim) || xdim < 2
        error('xdim must be an integer scalar >= 2.');
    end

    N_obs = numel(constr);
    Obs = struct('A', cell(1, N_obs), 'b', cell(1, N_obs));

    for i = 1:N_obs
        ci = constr{i};
        t = lower(string(ci.type));

        switch t
            case "pol"
                V = ci.vtList;
                if isempty(V)
                    Obs(i).A = zeros(xdim, 0);
                    Obs(i).b = zeros(0, 1);
                    continue;
                end
                if size(V,2) ~= 2
                    V = V.';   % allow 2xN input
                end

                x = V(:,1);
                y = V(:,2);

                % CCW test (shoelace sign)
                is_ccw = 0.5 * sum(x .* circshift(y, -1) - circshift(x, -1) .* y) > 0;

                N_vert = size(V,1);
                A2 = zeros(2, N_vert);
                b  = zeros(N_vert, 1);

                for j = 1:N_vert
                    p1 = V(j,:).';
                    p2 = V(mod(j, N_vert) + 1,:).';
                    dv = p2 - p1;

                    if is_ccw
                        n = [ dv(2); -dv(1) ];
                    else
                        n = [ -dv(2); dv(1) ];
                    end

                    nn = norm(n);
                    if nn > 0
                        n = n / nn;
                    end

                    A2(:,j) = n;
                    b(j)    = n' * p1;
                end

                Obs(i).A = padA(A2, xdim);
                Obs(i).b = b;

            case "cur"
                c = ci.center(:);
                r = ci.a;   % radius
                if numel(c) ~= 2
                    error('constr{%d}.center must be a 2-vector for type="cur".', i);
                end
                if ~isscalar(r) || r < 0
                    error('constr{%d}.a must be a nonnegative scalar radius for type="cur".', i);
                end

                m = 12;
                theta = (0:m-1) * (2*pi/m);

                A2 = [cos(theta); sin(theta)];  % 2 x m
                b  = A2.' * c + r;              % m x 1

                Obs(i).A = padA(A2, xdim);
                Obs(i).b = b;

            case "lin"
                a = ci.a(:);
                b = ci.b;

                if numel(a) == 2
                    A = zeros(xdim, 1);
                    A(1:2) = a;
                elseif numel(a) == xdim
                    A = a;
                else
                    error('constr{%d}.a must have length 2 or xdim=%d for type="lin".', i, xdim);
                end

                Obs(i).A = A;   % xdim x 1
                Obs(i).b = b;   % scalar

            otherwise
                error('Unsupported constraint type: %s (index %d).', t, i);
        end
    end
end

function A = padA(A2, xdim)
% Pad 2xm matrix A2 to xdim x m by adding zero rows.
    if isempty(A2)
        A = zeros(xdim, 0);
        return;
    end
    if size(A2,1) ~= 2
        error('padA expects a 2xm matrix.');
    end
    m = size(A2,2);
    A = zeros(xdim, m);
    A(1:2, :) = A2;
end
