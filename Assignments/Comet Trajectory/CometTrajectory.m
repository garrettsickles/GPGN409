function Comet_Trajectory(filename)
    [x, y, theta, s] = Import_Comet_Data(filename);
    Gx = [ones(length(x), 1) cosd(theta) sind(-theta)];
    Gy = [ones(length(y), 1) sind(theta) cosd(theta)];
    figure;
    Plot_Ellipse(LSWI(x,Gx,s,0), LSWI(y,Gy,s,0),...
        {'Color', [0 .76 1], 'LineWidth',1.5});
    Plot_Ellipse(LSWI(x,Gx,s,2.0),LSWI(y,Gy,s,2.0),...
        {'Color', [.08 .3 .5],'LineWidth',1.5});
    scatter(x,y,20,'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0 .5 .5],'LineWidth',1.0); hold on;
    for i=1:length(s)
        Plot_Circle([x(i), y(i), s(i)],{'k','LineWidth',1.5});
    end
    plot(0, 0,'-ko','LineWidth',1,'MarkerEdgeColor',[1 .5 0],...
        'MarkerFaceColor',[1 .5 0],'MarkerSize',10); hold on;
    Plot_Circle([0, 0, 1],...
        {'Color', [0 0.5 0],'LineWidth',2});
    Plot_Circle([0, 0, 5.20],...
        {'Color', [.9 .3 0],'LineWidth',2})
    title('Assignment 1: Inverting for the Trajectory of a Comet');
    xlabel('X (AU)'); ylabel('Y (AU)');
    axis([-7,13,-7,7]); axis equal;
    legend('Unweighted','Inverse Square','Observations','Uncertainty');
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
    W = diag(1./(w.^(n)));
    m = pinv((G')*(W')*(W)*(G))*(G')*(W')*(W)*(d);
end
function Plot_Ellipse(mx, my, param)
    t=0.0:0.001:2*3.1415926;
    plot(mx(1)+mx(2)*cos(t)-mx(3)*sin(t),...
        my(1)+my(2)*sin(t)+my(3)*cos(t),...
        param{:}); hold on;
end
function Plot_Circle(m, param)
    t=0.0:0.001:2*3.1415926;
    plot(m(1)+m(3).*cos(t),m(2)+m(3).*sin(t), param{:}); hold on;
end