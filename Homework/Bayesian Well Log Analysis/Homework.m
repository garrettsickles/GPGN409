function Homework()
    filename = 'data.txt';
    PlotPosterior(filename, 33);
    %LogPanel(filename);
end

function [depth,gamma,phi,rho,caliper] = Import_Well_Log(filename)
    data = importdata(filename);
    d = data.data;
    depth = d(:,1);
    gamma = d(:,2);
    phi = d(:,3);
    rho = d(:,4);
    caliper = d(:,5);
end

function [vp] = Vel_From_Phi(phi)
    vp = ((1-phi)./(6.64)+phi./(1.5)).^(-1);
end

function [phi] = Phi_From_Vel(vel)
    phi = (vel.^(-1)-(1/6.64))/(1/(1.5)-1./(6.64));
end

function [rho] = Rho_From_Phi(phi)
    rho = ((1-phi).*(2.71)+phi.*(1.0));
end

function [phi] = Phi_From_Rho(rho)
    phi = (rho-2.71)/(1.0-2.71);
end

function [vp] = Vel_From_Gamma(gamma)
    vp = 5.654 - 0.008.*gamma;
end

function [FF] = PlotPrior(filename, d)
    [depth,gamma,phi,rho,caliper] = Import_Well_Log(filename);
    
    figure;
    v = Vel_From_Gamma(gamma);
    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    [X,Y,Z] = meshgrid(x1,x2,x3);
    sigma = [std(v(:)) 0 0;0 std(rho(:)) 0; 0 0 std(phi(:))].*((abs(caliper(d)-5.999))/range(caliper(:)));
    F = mvnpdf([X(:) Y(:) Z(:)], [v(d) rho(d) phi(d)], sigma);
    F = reshape(F,length(x2),length(x1),length(x3));
	set(slice(X,Y,Z,F,[3.0:.5:8.0],[],[]),'edgecolor','none');
    colorbar;
    title(['Prior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    figure;
    FF = F;
    F(length(x2),:,:) = sum(FF(:,:,:),1);
    F(:,1,:) = sum(FF(:,:,:),2);
    F(:,:,1) = sum(FF(:,:,:),3);
    size(F)
	set(slice(X,Y,Z,F,[3.0],[3.0],[0.0]),'edgecolor','none');
    colorbar;
    title(['Prior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
end

function [FF] = PlotTheory()
    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    [X,Y,Z] = meshgrid(x1,x2,x3);
    [A,B] = meshgrid(x2,x3);
    F = zeros(length(x2),length(x1),length(x3));
    F = reshape(F,length(x2),length(x1),length(x3));
    sigma = [0.02 0; 0 0.01];
    for i=1:length(x1)
        phi = Phi_From_Vel(x1(i));
        F(:,i,:) = reshape(mvnpdf([A(:) B(:)], [Rho_From_Phi(phi) phi], sigma), length(x2), length(x3))';
    end
    figure;
    set(slice(X,Y,Z,F,[3.0:.5:8.0],[],[]),'edgecolor','none');
    colorbar;
    title(['Theorhetical Joint PDF']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    figure;
    FF = F;
    F(length(x2),:,:) = sum(FF(:,:,:),1);
    F(:,1,:) = sum(FF(:,:,:),2);
    F(:,:,1) = sum(FF(:,:,:),3);
    size(F)
	set(slice(X,Y,Z,F,[3.0],[3.0],[0.0]),'edgecolor','none');
    colorbar;
    title(['Theorhetical Joint PDF']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
end

function [FF] = PlotPosterior(filename, d)
    [depth,gamma,phi,rho,caliper] = Import_Well_Log(filename);
    
    F = PlotPrior(filename, d).*PlotTheory();
    
    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    [X,Y,Z] = meshgrid(x1,x2,x3);
    figure;
    set(slice(X,Y,Z,F,[3.0:.5:8.0],[],[]),'edgecolor','none');
    colorbar;
    title(['Posterior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    FF = F;
    F(length(x2),:,:) = sum(FF(:,:,:),1);
    F(:,1,:) = sum(FF(:,:,:),2);
    F(:,:,1) = sum(FF(:,:,:),3);
    figure;
    set(slice(X,Y,Z,F,[3.0],[3.0],[0.0]),'edgecolor','none');
    colorbar;
    title(['Posterior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    
    PRIOR = PlotPrior(filename, d);
    [m,ix] = max(PRIOR(:));
    [I1,I2,I3] = ind2sub(size(PRIOR),ix);
    
    POST = FF;
    [m,iix] = max(PRIOR);
    [II1,II2,II3] = ind2sub(size(POST),ix);
    length(PRIOR(1,:,3))
    figure;
    
    PRPDF = normpdf(x1,x1(I2),std(x1,PRIOR(I1,:,I3)));
    PSTPDF = normpdf(x1,x1(II2),std(x1,POST(II1,:,II3)));
    plot(x1,PRPDF, 'b'); hold on;
    plot(x1, PSTPDF, 'r');
    plot([x1(I2)-std(x1,PRIOR(I1,:,I3)) x1(I2)+std(x1,PRIOR(I1,:,I3))], [max(PRPDF) max(PRPDF)], 'b'); hold on;
    plot([x1(I2) x1(I2)], [0 max(PRPDF)], 'b'); hold on;
    plot([x1(II2)-std(x1,POST(II1,:,II3)) x1(II2)+std(x1,POST(II1,:,II3))], [max(PSTPDF) max(PSTPDF)], 'r'); hold on;
    plot([x1(II2) x1(II2)], [0 max(PSTPDF)], 'r'); hold on;
    xlabel('v (km/s)');
    legend('Prior','Posterior');
    title(['Prior vs. Posterior at ',num2str(depth(d)),' ft']);
end

function LogPanel(filename)
    [depth,gamma,phi,rho,caliper] = Import_Well_Log(filename);
    
    figure;
    subplot(1,6,1);
    plot(gamma(:),depth(:),'k');
    xlabel('GR (API)');
    ylabel('Depth (ft)')
    set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%d'));
    set(gca,'Ydir','reverse')
    subplot(1,6,2);
    plot(caliper(:),depth(:),'g');
    xlabel('CALI (in)');
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'Ydir','reverse')
    subplot(1,6,3);
    plot(phi(:),depth(:),'r');
    xlabel('\phi(.)');
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'Ydir','reverse')
    xlim([0.0 0.5]);
    subplot(1,6,4);
    plot(rho(:),depth(:),'b');
    xlabel('\rho(g/cc)');
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'Ydir','reverse');
    
    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    v = Vel_From_Gamma(gamma);
    
    vstd = std(v);
    Prior = zeros(length(depth),length(x1));
    for i=1:length(depth)
        Prior(i,:) = mvnpdf(x1(:),v(i),vstd);
    end
    subplot(1,6,5);
    surf(x1,depth,Prior,'edgecolor','none');
    title('Prior')
    xlabel('v (km/s)')
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'Ydir','reverse');
    axis tight;
    view(0,90);
    
    [X,Y,Z] = meshgrid(x1,x2,x3);
    [A,B] = meshgrid(x2,x3);
    THEORY = zeros(length(x2),length(x1),length(x3));
    THEORY = reshape(THEORY,length(x2),length(x1),length(x3));
    sigma1 = [0.02 0; 0 0.01];
    for i=1:length(x1)
        phi1 = Phi_From_Vel(x1(i));
        THEORY(:,i,:) = reshape(mvnpdf([A(:) B(:)], [Rho_From_Phi(phi1) phi1], sigma1), length(x2), length(x3))';
    end
    
    POST = zeros(length(depth),length(x1));
    size(POST)
    length(depth)
    sigma2 = [std(v(:)) 0 0;0 std(rho(:)) 0; 0 0 std(phi(:))];
    for i=1:length(depth)
        F = mvnpdf([X(:) Y(:) Z(:)], [v(i) rho(i) phi(i)], sigma2.*((abs(caliper(i)-5.9999))/range(caliper(:))));
        F = reshape(F,length(x2),length(x1),length(x3));
        F = F.*THEORY;
        [m,ix] = max(F(:));
        [I1,I2,I3] = ind2sub(size(F),ix);
        PDF = mvnpdf(x1(:),x1(I2),std(x1,F(I1,:,I3)));
        PDF = PDF./max(PDF(:));
        POST(i,:) = PDF;
    end
    
    subplot(1,6,6);
    surf(x1,depth,POST,'edgecolor','none');
    title('Posterior')
    xlabel('v (km/s)')
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'Ydir','reverse');
    axis tight;
    view(0,90);
end