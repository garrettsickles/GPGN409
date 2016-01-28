function HW1(filename)
    data = importdata(filename);
    figure;
    errorbar(data(:, 1),data(:, 2),data(:, 3),'rx');
    hold on;
    Weighted_Ballistic_Inversion(data(:, 1), data(:, 2), data(:, 3), 2.0);
    Weighted_Ballistic_Inversion(data(:, 1), data(:, 2), data(:, 3), 1.0);
    Weighted_Ballistic_Inversion(data(:, 1), data(:, 2), data(:, 3), 0.0);
    legend('Data/Uncertainty', 'r = 2.0', 'r = 1.0', 'r = 0.0');
    title('Ballistic Problem Inversion: Data and Multiple Inversions');
    xlabel('Horizontal Location (km)');
    ylabel('Vertical Location (km)');
end

function [m] = Weighted_Ballistic_Inversion(h, z, w, p)
    W = diag((w.^(-p)));
    G = ones(length(h), 3);
    G(:, 2) = G(:, 2) .* h;
    G(:, 3) = G(:, 3) .* (h .^ 2);
    m = pinv((G.')*(W.')*(W)*(G))*(G.')*(W.')*(W)*(z);
    r = roots([m(3) m(2) m(1)]);
    x = r(2):0.1:r(1);
    y = m(1).*ones(1, length(x)) + m(2).*x + m(3).*(x.^2);
    x0 = r(2);
    theta = atan(m(2)+(2.*m(3).*x0));
    vel = ((-9.81*(10^(-3)))./((2.0).*m(3).*(cos(theta).^2))).^0.5;
    display(['X initial (km): ', num2str(x0)]);
    display(['Theta (deg): ', num2str(theta.*180./(3.1415926))]);
    display(['Velocity (km/s): ', num2str(vel)]);
    plot(x, y);
    axis tight;
    hold on;
end