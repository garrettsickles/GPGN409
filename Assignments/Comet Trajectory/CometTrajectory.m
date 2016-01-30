% This function takes a filename containing comet trajectory data in the
% and plots the data as well as an inversion for the elliptical path
function Comet_Trajectory(filename)
    data = importdata(filename);
    x = data.data(:,1)
    y = data.data(:,2);
    angle = data.data(:,3);
    sigma = data.data(:,4);
    
    
    % Draw the data points
    figure;
    scatter(x, y, 'b', 'o')
    hold on;
    t = 0:0.001:(2*3.1415926);
    for i=1:length(sigma)
       plot(x(i)+sigma(i)*cos(t), y(i)+sigma(i)*sin(t), 'k');
       hold on;
    end
    axis equal;
    
end

% This function take a model and performs a generalized least squares
% inversion
% --- Input Parameters ---
% --- m: Model Vector
% --- Output Parameters ---
% --- d: Data Vector
% --- G: The Magical Matrix Operator
% --- W: Data Weighting Matrix
function [ m ] = Least_Squares_Inversion(d, G, W)
    m = pinv((G.')*(W.')*(W)*(G))*(G.')*(W.')*(W)*(d);
end

