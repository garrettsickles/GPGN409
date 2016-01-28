function HW1(filename)
    data = importdata(filename);
    m = Weighted_Ballistic_Inversion(data(:, 1), data(:, 2), data(:, 3));
    figure;
    errorbar(data(:, 1),data(:, 2),data(:, 3),'rx');
    hold on;
    r = roots([m(3) m(2) m(1)]);
    x = r(2):0.1:r(1);
    y = m(1).*ones(1, length(x)) + m(2).*x + m(3).*(x.^2);
    plot(x, y);
    h0 = r(2);
    theta = atan(m(2)+(2.*m(3).*h0));
    vel = ((-9.81*(10^(-3)))./((2.0).*m(3).*(cos(theta).^2))).^0.5;
    display(['X initial (Km): ', num2str(x0)]);
    display(['X initial (Km): ', num2str(h0)]);
    display(['Theta (deg): ', num2str(theta.*180./(3.1415926))]);
    display(['Velocity (Km/s): ', num2str(vel)]);
end

function [m] = Weighted_Ballistic_Inversion(h, z, w)
    W = diag((w.^(-1.0)));
    G = ones(length(h), 3);
    G(:, 2) = G(:, 2) .* h;
    G(:, 3) = G(:, 3) .* (h .^ 2);
    m = pinv((G.')*(W.')*(W)*(G))*(G.')*(W.')*(W)*(z);
end

function [m] = Unweighted_Ballistic_Inversion(h, z)
    Weighted_Ballistic_Inversion(h, z, ones(length(h), 1));
end