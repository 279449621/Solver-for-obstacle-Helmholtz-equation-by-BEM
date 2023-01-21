%反散射那本书96页远场的公式
function [uinf, nu, y, xhat, phi] = find_far_field(obj, sol, angles)
    N = obj.n;
    x   = obj.x;
    dx  = obj.dx;
    k   = obj.k;
    t = (0:2*N-1)'*pi/N;
    y = x(t);
    dy = dx(t);
    ds = sqrt(sum(dy.^2, 2));
    s = size(angles);
    if s(1) > s(2)
        angles = angles';
    end
    xhat = [cos(angles); sin(angles)];
    nu = unit_normal(obj, t);
    eta = sol.eta;
    phi = sol.phi;
    uinf = zeros(size(angles'));
    a = exp(-1i*pi/4) / sqrt(8*pi*k);
    for j = 1:length(xhat)
        xh = xhat(:,j);
        K = (k*nu*xh + eta) .* exp(-1i*k*y*xh) .* phi .* ds; 
        uinf(j) = a * sum(K * pi/N);
    end
end


function v = unit_normal(obj, t)
    dx = obj.dx;
    s = size(t);
    if s(2) > s(1)
        t = t';
    end
    v = dx(t) * [0 -1; 1 0];
    r = sqrt(sum(v.^2, 2));
    v = v ./ repmat(r, [1 2]);
end 