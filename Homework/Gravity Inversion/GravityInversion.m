function GravityInversion()
    % Get data from file
    filename = 'data.txt';
    data = importdata(filename);
    data = data.data;
    % Import the data
    x_data = data(:,1);
    z_data = data(:,2);
    dg_data = data(:,3);
    dgs_data = data(:,4);
    % Setup x space
    x_min = 0.0;
    x_max = 30.0;
    x_delta = 0.25;
    x = x_min:x_delta:x_max;
    % Setup z space
    z_min = 0.0;
    z_max = 3.0;
    z_delta = 0.25;
    z = z_min:z_delta:z_max;
    % Setup advised topography
    z_topo = (-1.0)+(x/40.0)+...
        (0.25*cos((3.1415926/2).*(x.^2)./(x_max-x_min)));
    % Useful lengths
    lxd = length(x_data);
    lxz = length(x)*length(z);
    % Setp meshgrid for plots and iterating
    [XZ, ZX] = meshgrid(x,z);
    % Convert meshgrid
    X = reshape(XZ',[1 lxz]);
    Z = reshape(ZX',[1 lxz]);
    
    % Initialize the K matrix
    K = zeros(lxd, lxz);
    for i=1:lxd
        for j=1:lxz
            K(i,j) = (6.67*10^(-2))*(10)*(2)*...
                (x_delta*z_delta)*(Z(j)-z_data(i))/...
                ((X(j)-x_data(i))^2+(Z(j)-z_data(i))^2)^(3/2);
        end
    end
    
    % Test model example
    tsm = zeros(length(x), length(z));
    tsm(46:50,4:8) = -100;
    
    figure;
    subplot(2,1,1);
    plot(x_data, GravityDataFromModel(K, reshape(tsm, lxz, 1)),...
        'color', 'b', 'LineWidth', 3);
    xlim([min(x),max(x)]);
    ylim([-400, 400]);
    title('G(m), \Delta \rho = -100 kg/m^{3}, y_{max} = 10 km');
    xlabel('x (km)');
    ylabel('\Delta g_{z} (km/s^{2})');
    
    subplot(2,1,2);
    plot(x, z_topo, 'color', 'b'); hold on;
    surf(XZ, ZX, reshape(tsm, length(x), length(z))',...
        'EdgeAlpha',0.0);
    xlim([min(x),max(x)]);
    ylim([min(z_topo).*1.2,max(z)]);
    title('m')
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    colormap winter;
    axis equal;
    axis tight;
    
    % Setup the inverse problem
    W_d = diag(1./(dgs_data));
    dg_model = zeros(lxz,1);
    iter = 30;
    S = logspace(-2,6,iter);
    
    W_mk = SparseLaplacian(lxz, length(x), length(z));

    [ dg, r_d, r_m ] = GetResiduals(x, z, dg_data, dg_model,...
        K, W_mk, W_d, S);
    
    index = OptFromResiduals(x_data, dg_data, dgs_data,...
        dg, K, S, r_d, r_m);
    
    PlotGravityModel(x_data,z_data,dg_data,dgs_data,...
        x,z,dg(:,1),K,GetCovariance(K,W_d,W_mk, S(1)),...
        z_topo, S(1));
    
    PlotGravityModel(x_data,z_data,dg_data,dgs_data,...
        x,z,dg(:,iter),K,GetCovariance(K,W_d,W_mk, S(iter)),...
        z_topo, S(iter));
    
    PlotGravityModel(x_data,z_data,dg_data,dgs_data,...
        x,z,dg(:,index),K,GetCovariance(K,W_d,W_mk, S(index)),...
        z_topo, S(index));
end

function [ C ] = GetCovariance(K, Wd, Wm, S)
    Wm = Wm./(S);
    C = pinv((K'*(Wd'*Wd)*K)+(Wm'*Wm));
end

function [ dg, rd, rm ] = GetResiduals(x, z, dgd, dgm, K, Wmk, Wd, S)
    n = length(S);
    dg = zeros(length(x)*length(z), n);
    rd = zeros(n, 1);
    rm = zeros(n, 1);
    for i=1:length(S)
        Wm = Wmk./(S(i));
        dg(:,i) = (pinv(K'*(Wd'*Wd)*K+Wm'*Wm)*...
            (K'*(Wd'*Wd)*dgd+Wm'*Wm*dgm))';
        rd(i) = (norm(Wd*(K*dg(:,i)-dgd)));
        rm(i) = (norm(Wmk*(dg(:,i)-dgm)));
        disp([num2str(i), ': ', num2str(S(i))]);
    end
end

function [ index ] = OptFromResiduals(x_data, dg_data, dgs_data,...
    dg, K, S, r_d, r_m)

    r_d = r_d./max(r_d(:));
    r_m = r_m./max(r_m(:));
    r_norm = (r_d.^2+r_m.^2).^(0.5);
    [junk, index] = min(r_norm);
    
    last = size(dg, 2);
    
    figure('Position',[400 400 600 400]);
    subplot(2,2,1);
    plot(r_m, r_d,...
        '.', 'MarkerFaceColor', 'k', 'color', 'k'); hold on;
    scatter(r_m(index), r_d(index),...
        'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); hold on;
    scatter(r_m(1), r_d(1),...
        'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    scatter(r_m(last), r_d(last),...
        'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    xlabel('||r_{m}||');
    ylabel('||r_{d}||');
    xlim([0.0 1.0]);
    ylim([0.0 1.0]);
    axis equal;
    title('Data and Model Residuals');
    
    subplot(2,2,2);
    plot(log10(S), sqrt(r_d.^2+r_m.^2),...
        '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
    plot(log10(S(index)), r_norm(index),...
        'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g'); hold on;
    plot(log10(S(1)), r_norm(1),...
        'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    plot(log10(S(last)), r_norm(last),...
        'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    xlabel('$\log_{10}(\sigma_{m})$','Interpreter','Latex');
    ylabel('$\sqrt{||r_{m}||^{2}+||r_{m}||^{2}}$','Interpreter','Latex');
    title('Minimized Residual');
    
    subplot(2,2,[3,4]);
    errorbar(x_data, dg_data, dgs_data,...
        'ro', 'color', 'k'); hold on;
    plot(x_data, GravityDataFromModel(K, dg(:,1)),...
        'color', 'r', 'LineWidth', 2); hold on;
    plot(x_data, GravityDataFromModel(K, dg(:,last)),...
        'color', 'b', 'LineWidth', 2); hold on;
    plot(x_data, GravityDataFromModel(K, dg(:,index)),...
        'color', 'g', 'LineWidth', 2); hold on;
    xlabel('X (km)');
    ylabel('\Delta g_{z} (km/s^{2})');
    axis tight;
    legend('Data Points','Min. Model','Max. Model','Opt. Model',...
        'Location','southeast');
    title('Raw Data with Certain Models Displayed');
end

function PlotGravityModel(x_data, z_data, dg_data, dgs_data,...
    x, z, dg, K, C, zs, sd)

    [XZ, ZX] = meshgrid(x,z);
    figure;
    subplot(3,1,1);
    errorbar(x_data, dg_data, dgs_data,...
        'ro', 'color', 'k'); hold on;
    plot(x_data, GravityDataFromModel(K, dg),...
        'color', 'b', 'LineWidth', 3);
    xlim([min(x),max(x)]);
    title(['$d, ', ' \sigma_{d} = ', num2str(sd), ', G(\tilde{m})$'],...
        'Interpreter','Latex');
    xlabel('x (km)');
    ylabel('\Delta g_{z} (km/s^{2})');
    
    subplot(3,1,2);
    plot(x, zs, 'color', 'b'); hold on;
    scatter(x_data, z_data, 'o', 'MarkerEdgeColor', 'k'); hold on;
    surf(XZ, ZX, reshape(dg, length(x), length(z))',...
        'EdgeAlpha',0.0);
    xlim([min(x),max(x)]);
    ylim([min(zs).*1.2,max(z)]);
    title('$\tilde{m}$','Interpreter','Latex');
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    
    subplot(3,1,3);
    plot(x, zs, 'color', 'b'); hold on;
    scatter(x_data, z_data, 'o', 'MarkerEdgeColor', 'k'); hold on;
    surf(XZ, ZX, reshape(diag(C), length(x), length(z))',...
        'EdgeAlpha',0.0);
    xlim([min(x),max(x)]);
    ylim([min(zs).*1.2,max(z)]);
    title('$diag(\tilde{C}_{m})$','Interpreter','Latex');
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
end

function [ d ] = GravityDataFromModel(G, m)
    d = G*m;
end

function [ L ] = SparseLaplacian(n, x, y)
    L = sparse(n,n);
    for i=1:n
        for j=1:n
            L(i,j) = LaplacianFilter2DAt(i,j,x,y);
        end
    end
end

function [ v ] = LaplacianFilter2DAt(r, c, x, y)
    if(r == c)
        v = -4;
    elseif(r-1 == c) % Below or to the left
        if(mod(r-1,x) == 0)
            v = 0;
        else
            v = 1;
        end
    elseif(c-1 == r) % Above or to the right
        if(mod(c-1,x) == 0)
            v = 0;
        else
            v = 1;
        end
    elseif(c-x == r || c+x == r)
        v = 1;
    else
        v = 0;
    end
end