function HorizonInterpolation()
    filename = 'data.txt';
    data = importdata(filename);
    data = data.data;
    x_data = data(:, 1);
    z_data = data(:, 2);
    s_data = data(:, 3);
    z_mean = mean(z_data);
    x_delta = 0.1;
    x_min = 0.0;
    x_max = 40.0;
    x = x_min:x_delta:x_max;
    z = ones(length(x), 1).*z_mean;
    K = zeros(length(x_data), length(x));
    for r=1:length(x_data)
        index = round((x_data(r)-x_min)./(x_delta));
        K(r, index) = 1;
    end
    z_data = z_data - z_mean;
    z = z - z_mean;
    S_m = 1.0;
    W_d = diag(1./(s_data.^2));
    W_m = diag(ones(length(x), 1)).*(-2)...
        + diag(ones(length(x)-1, 1), 1)...
        + diag(ones(length(x)-1, 1), -1);
    W_m = W_m./(x_delta*S_m);
    z = pinv(K'*(W_d'*W_d)*K+W_m'*W_m)*(K'*(W_d'*W_d)*z_data+W_m'*W_m*z);
    z = z + z_mean;
    z_data = z_data + z_mean;
    figure('Position', [0, 0, 800, 300]);
    errorbar(x_data, z_data, s_data, 'ro', 'color', 'k'); hold on;
    plot(x, z, 'color', 'b', 'LineWidth', 1);
    xlabel('X (km)');
    ylabel('Z (km)');
    axis tight;
    axis equal;
end