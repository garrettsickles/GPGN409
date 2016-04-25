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
    
    K = zeros(length(x_data), length(x_topo)); % Initialize the K matrix
    for r=1:length(x_data)
        index = round((x_data(r)-x_min)./(x_delta)); % Find closest point
        K(r, index) = 1; % Set entry to 1
    end
    
    size(LaplacianFilter2D(length(z_topo),length(x_topo)))
    
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

function [ B ] = LaplacianFilter2D(x,y)
    k = x*y;
    B = zeros(k);
    B = B + diag(ones(k,1).*(-4));
    B = B + diag(ones(k-1,1),1);
    B = B + diag(ones(k-1,1),-1);
    B = B + diag(ones(k-x,1),x);
    B = B + diag(ones(k-x,1),-x);
    for i=1:(y-1)
        j = i*x; 
        B(j,j+1) = 0;
        B(j+1,j) = 0;
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