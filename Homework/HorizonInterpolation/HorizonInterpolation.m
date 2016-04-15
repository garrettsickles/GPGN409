function HorizonInterpolation()
    filename = 'data.txt';
    data = importdata(filename);
    data = data.data;
    x_data = data(:, 1);
    z_data = data(:, 2);
    s_data = data(:, 3);

    x_delta = 0.1;
    x_min = 0.0;
    x_max = 40.0;
    x = x_min:x_delta:x_max;
    z = zeros(length(x), 1);
    K = zeros(length(x_data), length(x));
    for r=1:length(x_data)
        index = round((x_data(r)-x_min)./(x_delta));
        K(r, index) = 1;
    end

    S_m = 1.0;
    W_d = diag(1./(s_data.^2));
    W_mk = diag(ones(length(x), 1)).*(-2)...
        + diag(ones(length(x)-1, 1), 1)...
        + diag(ones(length(x)-1, 1), -1);
    W_m = W_mk./(x_delta*S_m);
    z = pinv(K'*(W_d'*W_d)*K+W_m'*W_m)*(K'*(W_d'*W_d)*z_data+W_m'*W_m*z);
    
    figure('Position', [0, 0, 800, 300]);
    errorbar(x_data, z_data, s_data, 'ro', 'color', 'k'); hold on;
    plot(x, z, 'color', 'b', 'LineWidth', 1);
    xlabel('X (km)');
    ylabel('Z (km)');
    axis tight;
    axis equal;
    
    
    iter = 100;
    r_d = zeros(iter, 1);
    r_m = zeros(iter, 1);
    S_mv = logspace(-3,3,iter);

    for i=1:iter
        S_m = S_mv(i);
        W_m = W_mk./(x_delta*S_m.^2);
        z = zeros(length(x), 1);
        z = pinv(K'*(W_d'*W_d)*K+W_m'*W_m)*(K'*(W_d'*W_d)*z_data+W_m'*W_m*z);
        r_d(i) = (norm(W_d*(K*z-z_data)));
        r_m(i) = (norm(W_mk*z));
    end

    figure;
    plot(r_d, r_m, 'o');
    
    figure;
    plot(S_mv, sqrt(r_d.^2+r_m.^2));
    
end