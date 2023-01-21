%提取参数
function [k, n, x, dx, ddx, f] = get_param(obj)
    k = obj.k;
    n = obj.n;
    x = obj.x;
    dx = obj.dx;
    ddx = obj.ddx;
    f = obj.f;
end