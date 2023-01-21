function [v, xv, yv] = find_u(x, y, obj, sol, t)
%计算Helmhotlz方程的解;
%x,y分别是考虑区域的恒坐标和纵坐标
    phi = sol.phi;
    n = length(t)/2;
    eta = obj.eta;
    [xv, yv] = make_param_vecs(x, y);
    [xvv, ~, s1] = make_param_vecs(xv, t);
    [yvv, tauv, ~] = make_param_vecs(yv, t);
    vn = [xvv, yvv];
    Koffbdy = (Loffbdy(vn, tauv, obj) + 1i*eta*Moffbdy(vn, tauv, obj))*pi/n;
    Koffbdy = reshape(Koffbdy, s1);
    v = Koffbdy * phi;
end

function v = Loffbdy(xv, tau, obj)
    [k, ~, x, dx, ~] = get_param(obj);
    xtau_t = x(tau) - xv;
    rotdxtau = dx(tau) * [0 -1; 1 0];
    nrmxt_tau =  sqrt(sum(xtau_t.^2, 2)); 
    v = -1i*k/4 * sum(rotdxtau.*xtau_t, 2) .* H1_1(k*nrmxt_tau) ./ nrmxt_tau;
end

function v = Moffbdy(xv, tau, obj)
    [k, ~, x, dx, ~] = get_param(obj);
    xt_tau = xv - x(tau);
    nrmxt_tau = sqrt(sum(xt_tau.^2, 2));
    dxtau = dx(tau);
    nrmdxtau = sqrt(sum(dxtau.^2, 2));
    v =  -1i/4 * H1_0(k*nrmxt_tau) .* nrmdxtau;
end
