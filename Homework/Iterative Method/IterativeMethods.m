function IterativeMethods()
    filename = 'data.txt';
    data = importdata(filename);
    x = data(:,1);
    z = data(:,2);
    
    d = [z(:)];
    G = [ones(length(x), 1) x(:) x(:).^2];
    
    a = -2:0.05:2;
    b = -2.0:0.05:2;
    c = -1.0:0.05:1.0;
    m_naught = [0 2.0 -0.5]';
    [mx, my, mz] = sphere;
    r = 0.06;
    
    SD = SteepestDescent(m_naught, d, G, 3);
    CG = ConjugateGradient(m_naught, d, G, 3);
    
    [J, X, Y, Z] = ObjectiveFunction(a, b, c, d, G);
    
    figure('Position',[400,400,600,450]);
    surf(mx*r+m_naught(1), my*2*r+m_naught(2), mz*r+m_naught(3),...
        'FaceColor',[1 0 0], 'EdgeAlpha', 0);
    hold on;
    set(slice(X,Y,Z,J,median(a),median(b),median(c)),'edgecolor','none');
    xlabel('m_{1} (km)');
    ylabel('m_{2} (-)');
    zlabel('m_{3} (km^{-1})');
    title('Objective Function Volume');
    colormap jet;
    axis tight;
    xlabel(colorbar, 'km');
    legend('m_{o}');
    view(-115, 17);
    
    figure('Position',[400,400,600,450]);
    set(slice(X,Y,Z,J,max(a),min(b),min(c)),'edgecolor','none');
    xlabel('m_{1} (km)');
    ylabel('m_{2} (-)');
    zlabel('m_{3} (km^{-1})');
    title('Objective Function Volume with Steepest Descent path');
    colormap jet;
    axis tight;
    xlabel(colorbar, 'km');
    hold on;
    plot3(SD(:,1),SD(:,2),SD(:,3),'r','LineWidth',7); hold on;
    hold on;
    surf(mx*r+m_naught(1), my*r+m_naught(2), mz*r+m_naught(3));
    axis equal;
    view(-115, 22);
    
    figure('Position',[400,400,600,450]);
    set(slice(X,Y,Z,J,max(a),min(b),min(c)),'edgecolor','none');
    xlabel('m_{1} (km)');
    ylabel('m_{2} (-)');
    zlabel('m_{3} (km^{-1})');
    title('Objective Function Volume with Conjugate Gradient path');
    colormap jet;
    axis tight;
    xlabel(colorbar, 'km');
    hold on;
    plot3(CG(:,1),CG(:,2),CG(:,3),'LineWidth',7); hold on;
    surf(mx*r+m_naught(1), my*r+m_naught(2), mz*r+m_naught(3));
    axis equal;
    view(-115, 22);
    
    figure('Position',[400,400,600,300]);
    BallisticPlot(SD, [x(:) z(:)]);
    xlabel('x (km)');
    ylabel('z (km)');
    title('Steepest Descent Method');
    legend('m^{o}','m^{1}','m^{2}','m^{3}','Data');
    
    figure('Position',[400,400,600,300]);
    BallisticPlot(CG, [x(:) z(:)]);
    xlabel('x (km)');
    ylabel('z (km)');
    title('Conjugate Gradient Method');
    legend('m^{o}','m^{1}','m^{2}','m^{3}','Data');
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
    for i=1:size(m,1)
        [x, y, p] = BallisticTrajectory(m(i,:));
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

function [x, y, p] = BallisticTrajectory(m)
    r = roots([m(3) m(2) m(1)]);
    x = min(r(:)):0.01:max(r(:));
    y = m(1).*ones(1, length(x)) + m(2).*x + m(3).*(x.^2);
    x0 = x(1);
    theta = atan(m(2)+(2.*m(3).*x0));
    v0 = ((-9.81*(10^(-3)))./((2.0).*m(3).*(cos(theta).^2))).^0.5;
    p = [x0 theta.*180./(3.1415926) v0];
end