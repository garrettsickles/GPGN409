function Comet_Trajectory(filename)
    % Import the data
    data = importdata(filename);
    x = data.data(:,1);
    y = data.data(:,2);
    angle = data.data(:,3);
    sigma = data.data(:,4);
    t = 0:0.001:(2*3.1415926);
    % Setup the Operator
    Gx = ones(length(x), 2);
    Gx(:,2) = cosd(angle);
    Gy = ones(length(y), 2);
    Gy(:,2) = sind(angle);
    % Unweighted Inversion
    mx = Least_Squares_Inversion(x, Gx, diag(sigma.^0));
    my = Least_Squares_Inversion(y, Gy, diag(sigma.^0));
    % Setup and run the inversion (^-1)
    mx1 = Least_Squares_Inversion(x, Gx, diag(sigma.^-1));
    my1 = Least_Squares_Inversion(y, Gy, diag(sigma.^-1));
    % Setup and run the inversion (^-2)
    mx2 = Least_Squares_Inversion(x, Gx, diag(sigma.^-2));
    my2 = Least_Squares_Inversion(y, Gy, diag(sigma.^-2));
    % Plot the trajectories
    figure;
    plot(mx(1)+mx(2)*cos(t), my(1)+my(2)*sin(t),...
        'Color', [0 .76 1], 'LineWidth',1.5);
    hold on;
    plot(mx1(1)+mx1(2)*cos(t), my1(1)+my1(2)*sin(t),...
        'Color', [.1 .5 .9], 'LineWidth',1.5);
    hold on;
    plot(mx2(1)+mx2(2)*cos(t), my2(1)+my2(2)*sin(t),...
        'Color', [.08 .3 .5],'LineWidth',1.5);
    hold on;
    % Plot the data points
    scatter(x,y,20,'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0 .5 .5],...
        'LineWidth',1.0);
    hold on;
    % Plot the error
    for i=1:length(sigma)
       plot(x(i)+sigma(i)*cos(t), y(i)+sigma(i)*sin(t),'k',...
           'LineWidth',2.0);
       hold on;
    end
    % Plot the sun
    plot(0, 0,'-ko',...
        'LineWidth',1,...
        'MarkerEdgeColor',[1 .5 0],...
        'MarkerFaceColor',[1 .5 0],...
        'MarkerSize',10);
    hold on;
    % Plot the Earth's Orbit
    plot(cos(t), sin(t),...
        'Color', [0 0.5 0],...
        'LineWidth',2);
    % Plot Jupiter's Orbit
    Jd = 5.20;
    plot(Jd*cos(t), Jd*sin(t),...
        'Color', [.9 .3 0],...
        'LineWidth',2);
    title('Assignment 1: Inverting for the Trajectory of a Comet');
    xlabel('X (AU)');
    ylabel('Y (AU)');
    axis([-7,13,-7,7])
    axis equal;
    legend('Unweighted', 'Inverse',...
        'Inverse Square', 'Observations',...
        'Uncertainty');
end

function [ m ] = Least_Squares_Inversion(d, G, W)
    m = pinv((G.')*(W.')*(W)*(G))*(G.')*(W.')*(W)*(d);
end

