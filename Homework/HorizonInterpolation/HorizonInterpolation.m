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
    W_d = diag(1./(s_data)); % Diagonalize the data uncertainties
    W_mk = diag(ones(length(x), 1)).*(-2)...
        + diag(ones(length(x)-1, 1), 1)...
        + diag(ones(length(x)-1, 1), -1); % Setup the interpolation matrix
    
    W_mk = W_mk./(x_delta.^2);
    
    iter = 101;
    z = zeros(length(x), iter);
    z_model = ones(length(x), 1);
    z_model(1:length(z_model)) = linspace(4.0,10.0,length(x));
    r_d = zeros(iter, 1);
    r_m = zeros(iter, 1);
    S_mv = logspace(-4,4,iter);
    
    % Uncomment the below section to add in a fault
%     x1 = 18;
%     x2 = 27;
%     ix1 = round((x1-x_min)./(x_delta));
%     ix2 = round((x2-x_min)./(x_delta));
%     W_mk(ix1, ix1-1) = 0;
%     W_mk(ix1-1, ix1) = 0;
%     W_mk(ix2, ix2+1) = 0;
%     W_mk(ix2+1, ix2) = 0;
%     z_model(1:(ix1-1)) = linspace(5,1,ix1-1);
%     z_model(ix1:ix2) = linspace(7.5,5,ix2-ix1+1);
%     z_model((ix2+1):length(z_model)) = linspace(4,10,length(z_model)-ix2);
    
    for i=1:iter
        W_m = W_mk./(S_mv(i));
        z(:,i) = pinv(K'*(W_d'*W_d)*K+W_m'*W_m)*(K'*(W_d'*W_d)*z_data+W_m'*W_m*z_model);
        r_d(i) = (norm(W_d*(K*z(:,i)-z_data)));
        r_m(i) = (norm(W_mk*(z(:,i)-z_model)));
    end
    
    r_d = r_d./max(r_d(:));
    r_m = r_m./max(r_m(:));
    r_norm = (r_d.^2+r_m.^2).^(0.5);
    [junk, index] = min(r_norm);

    W_m = W_mk./(S_mv(1));
    C_mt = pinv((K'*(W_d'*W_d)*K)+(W_m'*W_m));
    Rho = zeros(size(C_mt));
    for i=1:size(C_mt,1)
        for j=1:size(C_mt,2)
            Rho(i,j) = C_mt(i,j)/(C_mt(i,i).^(0.5).*C_mt(j,j).^(0.5));
        end
    end
    
    figure('Position',[400 400 600 400]);
    subplot(2,3,1);
    surf(C_mt,'EdgeAlpha',0.0);
    colormap jet;
    view([0 90])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'YDir','reverse');
    title('Min. Covariance');
    axis tight;
    axis equal;
    
    subplot(2,3,4);
    surf(Rho,'EdgeAlpha',0.0);
    colormap jet;
    view([0 90])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'YDir','reverse');
    title('Min. Correlation');
    axis tight;
    axis equal;
    
    W_m = W_mk./(S_mv(index));
    C_mt = pinv((K'*(W_d'*W_d)*K)+(W_m'*W_m));
    Rho = zeros(size(C_mt));
    for i=1:size(C_mt,1)
        for j=1:size(C_mt,2)
            Rho(i,j) = C_mt(i,j)/(C_mt(i,i).^(0.5).*C_mt(j,j).^(0.5));
        end
    end
    
    subplot(2,3,2);
    surf(C_mt,'EdgeAlpha',0.0);
    colormap jet;
    view([0 90])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'YDir','reverse');
    title('Opt. Covariance');
    axis tight;
    axis equal;
    
    subplot(2,3,5);
    surf(Rho,'EdgeAlpha',0.0);
    colormap jet;
    view([0 90])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'YDir','reverse');
    title('Opt. Correlation');
    axis tight;
    axis equal;
    
    W_m = W_mk./(S_mv(iter));
    C_mt = pinv((K'*(W_d'*W_d)*K)+(W_m'*W_m));
    Rho = zeros(size(C_mt));
    for i=1:size(C_mt,1)
        for j=1:size(C_mt,2)
            Rho(i,j) = C_mt(i,j)/(C_mt(i,i).^(0.5).*C_mt(j,j).^(0.5));
        end
    end
    
    subplot(2,3,3);
    surf(C_mt,'EdgeAlpha',0.0);
    colormap jet;
    view([0 90])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'YDir','reverse');
    title('Max. Covariance');
    axis tight;
    axis equal;
    
    subplot(2,3,6);
    surf(Rho,'EdgeAlpha',0.0);
    colormap jet;
    view([0 90])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'YDir','reverse');
    title('Max. Correlation');
    axis tight;
    axis equal;
    
    figure('Position',[400 400 600 400]);
    subplot(2,2,1);
    plot(r_m, r_d,...
        '.', 'MarkerFaceColor', 'k', 'color', 'k'); hold on;
    scatter(r_m(index), r_d(index),...
        'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); hold on;
    scatter(r_m(1), r_d(1),...
        'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    scatter(r_m(iter), r_d(iter),...
        'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    xlabel('||r_{m}||');
    ylabel('||r_{d}||');
    xlim([0.0 1.0]);
    ylim([0.0 1.0]);
    axis equal;
    title('Data and Model Residuals');
    
    subplot(2,2,2);
    plot(log10(S_mv), sqrt(r_d.^2+r_m.^2),...
        '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
    plot(log10(S_mv(index)), r_norm(index),...
        'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); hold on;
    plot(log10(S_mv(1)), r_norm(1),...
        'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    plot(log10(S_mv(iter)), r_norm(iter),...
        'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    xlabel('\sigma_{m}');
    ylabel('$\sqrt{||r_{m}||^{2}+||r_{m}||^{2}}$','Interpreter','Latex');
    title('Minimized Residual');
    
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
    legend('Data Points','Model','Min. Model','Max. Model','Opt. Model',...
        'Location','southeast');
    title('Raw Data with Certain Models Displayed');
    
    figure('Position',[400 400 600 400]);
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
    title('Model & Data Space');
    xlabel('$\log_{10}(\sigma_{m})$','Interpreter','Latex');
    ylabel('X (km)');
    zlabel('Z (km)');
    legend('Model Space','Min. Model','Max. Model','Opt. Model',...
        'Location','northwest');
end