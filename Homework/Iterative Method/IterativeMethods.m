function IterativeMethods()
    filename = 'data.txt';
    data = importdata(filename);
    x = data(:,1);
    z = data(:,2);
    
    d = [z(:)];
    G = [ones(length(x), 1) x(:) x(:).^2];
    
    a = -2:0.05:2;
    b = 0.5:0.05:2;
    c = -0.5:0.05:0.5;
    m0 = [0 1 -0.5]';
    [mx, my, mz] = sphere;
    r = 0.03;
    
    m1 = SteepestDescent(m0, d, G, 4);
    m2 = ConjugateGradient(m0, d, G, 3);
    
    ObjectiveFunction(a, b, c, d, G);
    plot3(m1(:,1),m1(:,2),m1(:,3),'r','LineWidth',7); hold on;
    surf(mx*r+m0(1), my*r+m0(2), mz*r+m0(3));
    axis equal;
    view(-103, 16);
    
    ObjectiveFunction(a, b, c, d, G);
    plot3(m2(:,1),m2(:,2),m2(:,3),'r','LineWidth',7); hold on;
    surf(mx*r+m0(1), my*r+m0(2), mz*r+m0(3));
    axis equal;
    view(-103, 16);
    
    BallisticPlot(m1, [x(:) z(:)]);
    BallisticPlot(m2, [x(:) z(:)]);
end

function [J, X, Y, Z] = ObjectiveFunction(x, y, z, d, G)
    [X, Y, Z] = meshgrid(x, y, z);
    OJ = @(m) ((m'*(G'*G)*m)./2-((G'*d)'*m)+(d'*d)./2);
    J = zeros(length(y), length(x), length(z));
    for i=1:length(x)
        for j=1:length(y)
            for k=1:length(z)
                J(j,i,k) = OJ([x(i) y(j) z(k)]');
            end
        end
    end
    
    figure;
    set(slice(X,Y,Z,J,max(x),min(y),min(z)),'edgecolor','none');
    colormap jet;
    axis tight;
    colorbar;
    hold on;
end

function [mv] = SteepestDescent(m, d, G, n)
    mv = zeros(n+1,size(G,2));
    mv(1,:) = m';
    for i=1:n
        A = G'*G;
        b = G'*d;
        m = mv(i,:)';
        p = b-A*m;
        mv(i+1,:) = (m+p*((p'*p)/(p'*A*p)))';
    end
end

function [mv] = ConjugateGradient(m, d, G, n)
    mv = zeros(n+1,size(G,2));
    mv(1,:) = m';
    A = G'*G;
    b = G'*d;
    r_old = b-A*m;
    p = r_old;
    for i=1:n
        m = mv(i,:)';
        a = (r_old'*r_old)/(p'*A*p);
        mv(i+1,:) = (m+a*p)';
        r_new = r_old - a*A*p;
        p = r_new + ((r_new'*r_new)/(r_old'*r_old))*p;
        r_old = r_new;
    end
end

function BallisticPlot(m, d)
    figure;
    
    for i=1:size(m,1)
        display(['Iteration #',num2str(i),':'])
        [x, y] = BallisticTrajectory(m(i,:));
        plot(x, y,...
            '--','LineWidth', 1.5);
        hold on;
    end
    scatter(d(:,1), d(:,2),...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',[0 .5 .5],...
        'LineWidth',1.0);
    axis tight;
    hold on;
end

function [x, y] = BallisticTrajectory(m)
    r = roots([m(3) m(2) m(1)]);
    x = r(2):0.1:r(1);
    y = m(1).*ones(1, length(x)) + m(2).*x + m(3).*(x.^2);
    x0 = r(2);
    theta = atan(m(2)+(2.*m(3).*x0));
    v0 = ((-9.81*(10^(-3)))./((2.0).*m(3).*(cos(theta).^2))).^0.5;
    display(['X initial (km): ', num2str(x0)]);
    display(['Theta (deg): ', num2str(theta.*180./(3.1415926))]);
    display(['Velocity (km/s): ', num2str(v0)]);
end