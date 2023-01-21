function [sol, obj] = solve_for_phi_new_bdy(sol, obj, objnew)
%求解不同边界下对应的解phi，用于求解反问题
    solver = sol.solver;
    t = solver.t;
    g = getg(t, objnew);
    solver.g = g;
    phi = solver.OP \ solver.g;
    % Update solution structure
    sol.solver = solver;
    sol.phi = phi;
    sol.phiint = @(t) phiint(t, obj, sol);
    sol.find_u = @(x, y, t) find_u(x, y, obj, sol, t);
    sol.find_u_list = @(x, t) find_u_list(x, obj, sol, t);
    sol.find_uinf = @(angles) find_far_field(obj, sol, angles);
    sol.update_bdy = @(objnew) solve_for_phi_new_bdy(sol, obj, objnew);
    obj = objnew;
end
function g = getg(t, obj)
    s = size(t);
    if s(2) > s(1)
        t = t';
    end
    [~,~,~,~,~,f] = get_param(obj);
    g = 2*f(t);
end

