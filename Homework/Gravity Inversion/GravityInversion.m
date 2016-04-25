function GravityInversion()
    filename = 'data.txt';
    data = importdata(filename);
    data = data.data;
    
    x_data = data(:,1);
    z_data = data(:,2);
    dg_data = data(:,3);
    dgs_data = data(:,4);
    
    x_min = 0.0;
    x_max = 30.0;
    x_delta = 0.2;
    x_topo = x_min:x_delta:x_max;
    z_topo = (-1.0)+(x_topo/40.0)+(0.25*cos((3.1415926/2).*(x_topo.^2)./(x_max-x_min)));
    
    figure;
    subplot(2,1,1);
    errorbar(x_data, dg_data, dgs_data,...
        'ro', 'color', 'k');
    xlim([x_min,x_max]);
    
    subplot(2,1,2);
    plot(x_topo, z_topo, 'color', 'b'); hold on;
    scatter(x_data, z_data, 'o', 'MarkerEdgeColor', 'k');
    xlim([x_min,x_max]);
end

