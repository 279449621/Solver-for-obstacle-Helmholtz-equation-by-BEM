function sol = dos_solve(obj)
    [phi, solver] = solve_for_phi(obj);
    sol.phi = phi;
    sol.solver = solver;
    sol.eta = obj.eta;
end
%反散射那本书92页公式L_2(t,t)=L(t,t)
function [v,si] = Ldiag(t,tau,obj)
    [~,~,~,dx,ddx] = get_param(obj);
    si = (t==tau);
    tsi = t(si);
    dxt = dx(tsi);
    rotdxt = ddx(tsi)*[0,-1;1,0];
    nrmdxt2 = sum(dxt.^2,2);
    v = 1/(2*pi)*sum(dxt.*rotdxt,2)./nrmdxt2;
end
%反散射那本书92页公式L(t,tau)
function v = L(t,tau,obj)
    [k, ~, x, dx, ~] = get_param(obj);
    xtau_t = x(tau)-x(t);
    rotdxtau = dx(tau)*[0,-1;1,0];
    nrmxt_tau = sqrt(sum(xtau_t.^2,2));
    v = 1i*k/2*sum(rotdxtau.*xtau_t,2).*H1_1(k*nrmxt_tau)./nrmxt_tau;
    [diagv,si] = Ldiag(t,tau,obj);
    v(si) = diagv;
end
%反散射那本书92页公式L1(t,tau)
function v = L1(t,tau,obj)
    [k, ~, x, dx, ~] = get_param(obj);
    xt_tau = x(t)-x(tau);
    rotdxtau = dx(tau)*[0,-1;1,0];
    nrmxt_tau =  sqrt(sum(xt_tau.^2, 2)); 
    v = k/(2*pi) * sum(rotdxtau.*xt_tau, 2) .* J_1(k*nrmxt_tau) ./ nrmxt_tau;
    [diagv, si] = Ldiag(t, tau, obj);
    v(si) = zeros(size(diagv));
end
%反散射那本书92页公式L2(t,tau)
function v = L2(t, tau, obj)
    v = L(t, tau, obj) - L1(t, tau, obj).*log(4*sin((t-tau)/2).^2);
    [diagv, si] = Ldiag(t, tau, obj);
    v(si) = diagv;
end
%反散射那本书92页公式M2(t,t)
function [v, si] = Mdiag(t, tau, obj)
    [k,~,~,dx,~] = get_param(obj);
    si = (t == tau);
    tsi = t(si);
    dxt = dx(tsi);
    nrmdxt  = sqrt(sum(dxt.^2, 2));
    C = 0.57721;
    v = (1i/2 - C/pi - 1/pi * log(k/2 * nrmdxt)).*nrmdxt;
end
%反散射那本书92页公式M(t,tau)
function v = M(t, tau, obj)
    [k,~, x, dx,~] = get_param(obj);
    xt_tau = x(t) - x(tau);
    nrmxt_tau = sqrt(sum(xt_tau.^2,2));
    dxtau = dx(tau);
    nrmdxtau = sqrt(sum(dxtau.^2,2));
    v = 1i/2*H1_0(k*nrmxt_tau).* nrmdxtau;
end
%反散射那本书92页公式M1(t,tau)
function v = M1(t, tau, obj)
    [k,~,x,dx,~] = get_param(obj);
    xt_tau = x(t) - x(tau);
    nrmxt_tau = sqrt(sum(xt_tau.^2, 2));
    dxtau = dx(tau);
    nrmdxtau = sqrt(sum(dxtau.^2, 2));
    v = -1/(2*pi) * J_0(k*nrmxt_tau) .* nrmdxtau;
end
%反散射那本书92页公式M2(t,tau)
function v = M2(t, tau, obj)
    v = M(t, tau, obj) - M1(t, tau, obj) .* log(4*sin((t-tau) / 2).^2);
    [diagv, si] = Mdiag(t, tau, obj);
    v(si) = diagv;
end
%反散射那本书93页公式K(t,tau)
function [K1, K2] = getK(t, obj)
    [~,n,~,~,~] = get_param(obj);
    eta = obj.eta;
    tau = (0:(2*n-1)) * pi/n;
    [t, tau, s] = make_param_vecs(t, tau);
    K1 = L1(t, tau, obj) + 1i*eta*M1(t, tau, obj);
    K2 = L2(t, tau, obj) + 1i*eta*M2(t, tau, obj);
    K1 = reshape(K1, s);
    K2 = reshape(K2, s);
end
%反散射那本书93页公式g(t)
function g = getg(t, obj)
    s = size(t);
    if s(2) > s(1)
        t = t';
    end
    [~,~,~,~,~,f] = get_param(obj);
    g = 2*f(t);
end
%设定矩阵方程组的系数
function solver = setup_solver(obj)
    [~,n,~,~,~] = get_param(obj);
    t = (0:(2*n-1))' * pi/n;
    [K1, K2] = getK(t, obj);
    R = dos_quad(n, t);
    A = R.*K1 + pi/n*K2;
    g = getg(t, obj);
    OP = eye(length(A)) - A;
    solver.OP = OP;
    solver.t = t;
    solver.g = g;
    solver.R = R;
    solver.A = A;
end
%求解线性方程组，得到反散射那本书93页结果Phi
function [phi, solver] = solve_for_phi(obj)
    solver = setup_solver(obj);
    phi = solver.OP \ solver.g; 
end