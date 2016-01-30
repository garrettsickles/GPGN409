% This function takes a filename containing comet trajectory data in the
% and plots the data as well as an inversion for the elliptical path
function Comet_Trajectory(filename)
    data = importdata(filename);
    x = data.data(:,1)
    y = data.data(:,2);
    angle = data.data(:,3);
    sigma = data.data(:,4);
    Gx = ones(length(x), 2);
    Gx(:,2) = cosd(angle);
    mx = Least_Squares_Inversion(x, Gx, diag(sigma.^-2))
    Gy = ones(length(y), 2);
    Gy(:,2) = sind(angle);
    my = Least_Squares_Inversion(y, Gy, diag(sigma.^-2))
    figure;
    scatter(x, y, 'b', 'o')
    hold on;
    t = 0:0.001:(2*3.1415926);
    for i=1:length(sigma)
       plot(x(i)+sigma(i)*cos(t), y(i)+sigma(i)*sin(t), 'k');
       hold on;
    end
    plot(mx(1)+mx(2)*cos(t), my(1)+my(2)*sin(t), 'b');
    axis equal;
    
end

function [ m ] = Least_Squares_Inversion(d, G, W)
    m = pinv((G.')*(W.')*(W)*(G))*(G.')*(W.')*(W)*(d);
end

