%将t和tau变为对应长度的向量
function [t, tau, matsize] = make_param_vecs(t, tau)
    %t为边界离散化，tau为积分离散化
    s = size(t);
    %t变为列向量
    if (s(2) > s(1))
        t = t';
    end
    s = size(tau);
    %tau变为行向量
    if (s(1) > s(2))
        tau = tau';
    end
    s = size(t);
    t = repmat(t, size(tau));
    tau = repmat(tau, s);
    matsize = size(t);
    vecsize = [matsize(1)*matsize(2),1];
    tau = reshape(tau,vecsize);
    t = reshape(t,vecsize);
end