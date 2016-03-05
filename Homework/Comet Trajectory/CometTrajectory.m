function Comet_Trajectory(filename)
    [x, y, theta, s] = Import_Comet_Data(filename);
    G = [ones(length(s), 1) cosd(theta) sind(theta)];
    figure('Position', [100, 100, 700, 500]);
    Trajectory(LSWI(x,G,s,0), LSWI(y,G,s,0),strcat('0-',filename),...
        {'Color', [0 .76 1], 'LineWidth',1.5});
    Trajectory(LSWI(x,G,s,2.0),LSWI(y,G,s,2.0),strcat('2-',filename),...
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
    title(['Assignment 1: Inverting for the Trajectory of a Comet in ',...
        filename]);
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
function [ h, k, a, b, p ] = Ellipse_Parameters(mx, my, filename)
    file = fopen(['Solved-', filename], 'w');
    h = mx(1); k = my(1);
    p = atan2d(my(2), mx(2));
    a = mx(2)./cosd(p); b = my(3)./cosd(p);
    fprintf(file, 'mx: %f, %f, %f\nmy: %f, %f, %f\n',...
        mx(1), mx(2), mx(3), my(1), my(2), my(3));
    fprintf(file,'h: %f\nk: %f\na: %f\nb: %f\nPhi: %f\n',...
        h, k, a, b, p);
    fclose(file);
end
function Trajectory(mx, my, id, param)
    t=0.0:0.001:2*3.1415926;
    Ellipse_Parameters(mx, my, id);
    plot(mx(1)+mx(2)*cos(t)+mx(3)*sin(t),...
        my(1)+my(2)*cos(t)+my(3)*sin(t),...
        param{:}); hold on;
end
function Plot_Circle(m, param)
    t=0.0:0.001:2*3.1415926;
    plot(m(1)+m(3).*cos(t),m(2)+m(3).*sin(t), param{:}); hold on;
end