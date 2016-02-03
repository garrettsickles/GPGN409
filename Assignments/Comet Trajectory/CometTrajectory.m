% x = h + a * cos(t) * cos(p) - b * sin(t) * sin(p)
% y = k + b * sin(t) * cos(p) + a * cos(t) * sin(p)

function Comet_Trajectory(filename)
    % Import the data
    [x, y, theta, s] = Import_Comet_Data(filename);
    % Setup the Operator
    Gx = [ones(length(x), 1) cosd(theta) sind(-theta)];
    Gy = [ones(length(y), 1) sind(theta) cosd(theta)];
    % Plot the trajectories
    figure;
    Plot_Ellipse(LSI(x,Gx), LSI(y,Gy),...
        {'Color', [0 .76 1], 'LineWidth',1.5});
    Plot_Ellipse(LSWI(x,Gx,s,-1), LSWI(y,Gy,s,-1),...
        {'Color', [.1 .5 .9], 'LineWidth',1.5});
    Plot_Ellipse(LSWI(x,Gx,s,-2),LSWI(y,Gy,s,-2),...
        {'Color', [.08 .3 .5],'LineWidth',1.5});
    % Plot the data points
    scatter(x,y,20,'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0 .5 .5],...
        'LineWidth',1.0);
    hold on;
    % Plot the error
    for i=1:length(s)
        Plot_Circle([x(i), y(i), s(i)],...
            {'k','LineWidth',1.5});
    end
    % Plot the sun
    plot(0, 0,'-ko',...
        'LineWidth',1,...
        'MarkerEdgeColor',[1 .5 0],...
        'MarkerFaceColor',[1 .5 0],...
        'MarkerSize',10);
    hold on;
    % Plot the Earth's Orbit
    Plot_Circle([0, 0, 1],...
        {'Color', [0 0.5 0],'LineWidth',2});
    % Plot Jupiter's Orbit
    Plot_Circle([0, 0, 5.20],...
        {'Color', [.9 .3 0],'LineWidth',2})
    % Setup Plot Stuff
    title('Assignment 1: Inverting for the Trajectory of a Comet');
    xlabel('X (AU)');
    ylabel('Y (AU)');
    axis([-7,13,-7,7])
    axis equal;
    legend('Unweighted', 'Inverse',...
        'Inverse Square', 'Observations',...
        'Uncertainty');
end

function [ x, y, theta, sigma ] = Import_Comet_Data(filename)
    data = importdata(filename);
    d = data.data;
    x = d(:,1);
    y = d(:,2);
    theta = d(:,3);
    sigma = d(:,4);
end

function [ m ] = LSWI(d, G, w, n)
    W = diag(w.^(n));
    m = pinv((G.')*(W.')*(W)*(G), 0.0001)*(G.')*(W.')*(W)*(d);
end

function [ m ] = LSI(d, G)
    m = pinv((G.')*(G))*(G.')*(d);
end

function [ h, k, a, b, p ] = Ellipse_From_Coefficients(mx, my)
    h = mx(1);
    k = my(1);
    p = atan2d(-mx(3), my(2));
    a = mx(2)./cosd(p);
    b = my(2)./cosd(p);
end

function Plot_Ellipse(mx, my, param)
    t=0.0:0.001:2*3.1415926;
    plot(mx(1)+mx(2)*cos(t)-mx(3)*sin(t),...
        my(1)+my(2)*sin(t)+my(3)*cos(t),...
        param{:});
    hold on;
end

function Plot_Circle(m, param)
    t=0.0:0.001:2*3.1415926;
    plot(m(1)+m(3).*cos(t),m(2)+m(3).*sin(t), param{:});
    hold on;
end
