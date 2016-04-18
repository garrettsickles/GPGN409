function HorizonInterpolation()
    % Import the data
    filename = 'data.txt';
    data = importdata(filename);
    data = data.data;
    x_data = data(:, 1); % X data
    z_data = data(:, 2); % Z data
    s_data = data(:, 3); % Std. Dev. in each Z measurement
    
    % Setup the K matrix
    x_delta = 0.1; % Distance between discrete X measurements
    x_min = 0.0; % Min X
    x_max = 40.0; % Max X
    x = x_min:x_delta:x_max; % Create X vector
    K = zeros(length(x_data), length(x)); % Initialize the K matrix
    for r=1:length(x_data)
        index = round((x_data(r)-x_min)./(x_delta)); % Find closest point
        K(r, index) = 1; % Set entry to 1
    end
    
    % Setup the inverse problem
    W_d = diag(1./(s_data.^2)); % Diagonalize the data uncertainties
    W_mk = diag(ones(length(x), 1)).*(-2)...
        + diag(ones(length(x)-1, 1), 1)...
        + diag(ones(length(x)-1, 1), -1); % Setup the interpolation matrix
    
    iter = 101;
    z = zeros(length(x), iter);
    z_model = ones(length(x), 1);
    r_d = zeros(iter, 1);
    r_m = zeros(iter, 1);
    S_mv = logspace(-2,2,iter);

    for i=1:iter
        W_m = W_mk./(S_mv(i));
        z(:,i) = pinv(K'*(W_d'*W_d)*K+W_m'*W_m)*(K'*(W_d'*W_d)*z_data+W_m'*W_m*z_model);
        r_d(i) = (norm(W_d*(K*z(:,i)-z_data)));
        r_m(i) = (norm(W_mk*(z(:,i)-z_model)));
    end
    
    r_norm = (r_d.^2+r_m.^2).^(0.5);
    [junk, index] = min(r_norm);
    
    figure;
    subplot(2,2,1);
    plot(r_m, r_d,...
        'o', 'MarkerFaceColor', 'k', 'color', 'k'); hold on;
    scatter(r_m(index), r_d(index),...
        'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); hold on;
    scatter(r_m(1), r_d(1),...
        'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    scatter(r_m(iter), r_d(iter),...
        'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    axis equal;
    xlabel('||r_{m}||');
    ylabel('||r_{d}||');
    title('Data and Model Residual Comparison');
    
    subplot(2,2,2);
    plot(log10(S_mv), sqrt(r_d.^2+r_m.^2),...
        'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
    plot(log10(S_mv(index)), r_norm(index),...
        'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); hold on;
    plot(log10(S_mv(1)), r_norm(1),...
        'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    plot(log10(S_mv(iter)), r_norm(iter),...
        'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    xlabel('\sigma_{m}');
    ylabel('$\sqrt{||r_{m}||^{2}+||r_{m}||^{2}}$','Interpreter','Latex');
    title('Minimized Residual with respect to Uncertainty');
    
    subplot(2,2,[3,4]);
    errorbar(x_data, z_data, s_data, 'ro', 'color', 'k'); hold on;
    scatter(x, z_model, '.',...
        'MarkerEdgeColor','m',...
        'MarkerFaceColor','m'); hold on;
    plot(x, z(:,1), 'color', 'r', 'LineWidth', 2); hold on;
    plot(x, z(:,iter), 'color', 'b', 'LineWidth', 2); hold on;
    plot(x, z(:,index), 'color', 'g', 'LineWidth', 2); hold on;
    xlabel('X (km)');
    ylabel('Z (km)');
    axis tight;
    axis equal;
    legend('Data Points','Model','Min. Model','Max. Model','Opt. Model');
    title('Raw Data with Certain Models Displayed');
    
    figure;
    [X, Y] = meshgrid(log10(S_mv), x);
    surf(X, Y, z, 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.1);
    colormap('gray');
    hold on;
    plot3(log10(S_mv(1)).*ones(length(x),1),x,z(:,1),...
        'Color', 'r', 'LineWidth', 4); hold on;
    plot3(log10(S_mv(iter)).*ones(length(x),1),x,z(:,iter),...
        'Color', 'b', 'LineWidth', 4); hold on;
    plot3(log10(S_mv(index)).*ones(length(x),1),x,z(:,index),...
        'Color', 'g', 'LineWidth', 4); hold on;
    view([43, 54]);
    title('Model Space with Specific Models');
    xlabel('X (km)');
    ylabel('Uncertainty');
    zlabel('Z (km)');
    legend('Model Space','Min. Model','Max. Model','Opt. Model');
end