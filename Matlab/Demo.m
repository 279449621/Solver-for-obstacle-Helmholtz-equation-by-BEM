clear;clc;
%定义边界曲线的方程(2D情形)
x   = @(t) [-0.65 + cos(t) + 0.65*cos(2*t),  1.5*sin(t)];
dx  = @(t) [      - sin(t) - 1.30*sin(2*t),  1.5*cos(t)];
ddx = @(t) [      - cos(t) - 2.60*cos(2*t), -1.5*sin(t)];
%考虑方程Delta u+k^2*u=0
k = 1;
%定义入射波
inc_ang = 0;
inc_dir = [cos(pi * inc_ang/180); sin(pi * inc_ang/180)];
%这里为边界条件，Dirichlet边界条件下，为负入射
f = @(t) -exp(1i*k*x(t)*inc_dir);
%定义离散点数目
n = 32;%积分分划方法: t = (0:(2*n-1))' * pi/n;
%定义单双势的权重
eta = k;
%定义传入变量:
obj.x = x;
obj.dx = dx;
obj.ddx = ddx;
obj.f = f;
obj.k = k;
obj.n = n;
obj.eta = eta;
%求解离散后矩阵方程的A x b满足Ax=b;
sol = dos_solve(obj);
phi = sol.phi;
%求解全局解
xn = 200;
yn = 200;
xs = linspace(-8,8,xn);
ys = linspace(-8,8,yn);
ts = sol.solver.t;
[v, xv, yv] = find_u(xs, ys,obj,sol,ts);
%求解远场
angles = linspace(0, 2*pi, 129)';
angles = angles(1:end-1);
[uinf, nu, y, xhat, phi] = find_far_field(obj, sol, angles);
%可视化
numperiods = 2;
c0 = 1;
w = k*c0;
times = linspace(0, numperiods*(2*pi/w), 16)';
timesnew = linspace(times(1), times(end), 128);
colormap(gray);
vm = reshape(v, xn, yn)';
ui = exp(1i*k*[xv, yv]*inc_dir);
ui = reshape(ui, xn, yn)';
for j = 1:length(timesnew)
    subplot(1,2,1);
    tds = exp(1i*w*timesnew(j))*vm;
    tdi = exp(1i*w*timesnew(j))*ui;
    imagesc(xs, ys, real(tds + tdi), [-2 2]);
    title('Total');
    subplot(1,2,2);
    imagesc(xs, ys, real(tds), [-1.5 1.5]);
    title('Scattered');
    pause(0.1);
end