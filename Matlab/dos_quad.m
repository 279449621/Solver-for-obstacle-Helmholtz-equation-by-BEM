function R = dos_quad(n, tn)
%反散射93页R_j^(n)公式
    s = size(tn);
    if s(2) > s(1)
        tn = tn';
    end
    tj = (0:(2*n-1))*pi/n;
    m = (1:(n-1))';
    R = zeros(length(tn), length(tj));
    for i = 1:length(tn)
        ti = tn(i);
        M = cos(m*(ti - tj)) ./ repmat(m, size(tj));
        R(i,:) = -2*pi/n*sum(M, 1) - pi/n^2 * cos(n*(ti - tj));
    end
    
end