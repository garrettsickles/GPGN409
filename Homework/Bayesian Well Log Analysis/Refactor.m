function Refactor()
    filename = 'data.txt';
    d = 3000;
    
    Prior = PlotPrior(filename, d);
    Theory = PlotTheory();
    PlotPosterior(filename, Prior, Theory, d);
    PlotPosterior(filename, Prior.*Theory, Theory, d);
    LogPanel(filename);
end

function [depth,gamma,phi,rho,caliper] = ImportWellLog(filename)
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

function [X, Y, Z, Prior, PriorP] = CalculatePrior(vel, rho, phi, v, r, p, sigma, d)
    [X,Y,Z] = meshgrid(v,r,p);
    F = mvnpdf([X(:) Y(:) Z(:)], [vel(d) rho(d) phi(d)], sigma);
    Prior = reshape(F,length(r),length(v),length(p));
    PriorP = Prior;
    PriorP(length(r),:,:) = sum(Prior(:,:,:),1);
    PriorP(1,:,:) = PriorP(length(r),:,:);
    PriorP(:,length(v),:) = sum(Prior(:,:,:),2);
    PriorP(:,1,:) = PriorP(:,length(v),:);
    PriorP(:,:,length(p)) = sum(Prior(:,:,:),3);
    PriorP(:,:,1) = PriorP(:,:,length(p));
end

function [F] = PlotPrior(filename, d)
    [depth,gamma,phi,rho,caliper] = ImportWellLog(filename);
    vel = Vel_From_Gamma(gamma);
    
    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    sigma = [std(vel(:)) 0 0;0 std(rho(:)) 0.0; 0 0.0 std(phi(:))].*10.*((abs(caliper(d)-5.999))/range(caliper(:)));
    [X, Y, Z, F, FF] = CalculatePrior(vel(:), rho(:), phi(:), x1(:), x2(:), x3(:), sigma, d);
	
    figure;
    set(slice(X,Y,Z,F,3.0:.5:8.0,[],[]),'edgecolor','none');
    colorbar;
    title(['Prior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    
    figure;
	set(slice(X,Y,Z,FF,3.0,3.0,0.0),'edgecolor','none');
    colorbar;
    title(['Prior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
end

function [X, Y, Z, Th, ThP] = CalculateTheory(v, r, p, sigma)
    [X,Y,Z] = meshgrid(v,r,p);
    [A,B] = meshgrid(r,p);
    Th = zeros(length(r),length(v),length(p));
    Th = reshape(Th,length(r),length(v),length(p));
    for i=1:length(v)
        phi = Phi_From_Vel(v(i));
        Th(:,i,:) = reshape(mvnpdf([A(:) B(:)], [Rho_From_Phi(phi) phi], sigma), length(r), length(p))';
    end
    ThP = Th;
    ThP(length(r),:,:) = sum(Th(:,:,:),1);
    ThP(1,:,:) = ThP(length(r),:,:);
    ThP(:,length(v),:) = sum(Th(:,:,:),2);
    ThP(:,1,:) = ThP(:,length(v),:);
    ThP(:,:,length(p)) = sum(Th(:,:,:),3);
    ThP(:,:,1) = ThP(:,:,length(p));
end

function [F] = PlotTheory()
    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    [X, Y, Z, F, FF] = CalculateTheory(x1(:), x2(:), x3(:), [0.02 0; 0 0.01]);
    
    figure;
    set(slice(X,Y,Z,F,3.0:.5:8.0,[],[]),'edgecolor','none');
    colorbar;
    title('Theorhetical Joint PDF');
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    
    figure;
	set(slice(X,Y,Z,FF,3.0,3.0,0.0),'edgecolor','none');
    colorbar;
    title('Theorhetical Joint PDF');
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
end

function [X, Y, Z, Post, PostP] = CalculatePosterior(Prior, Theory, v, r, p)
    [X,Y,Z] = meshgrid(v,r,p);
    Post = Prior.*Theory;
    
    PostP = Post;
    PostP(length(r),:,:) = sum(Post(:,:,:),1);
    PostP(1,:,:) = PostP(length(r),:,:);
    PostP(:,length(v),:) = sum(Post(:,:,:),2);
    PostP(:,1,:) = PostP(:,length(v),:);
    PostP(:,:,length(p)) = sum(Post(:,:,:),3);
    PostP(:,:,1) = PostP(:,:,length(p));
end

function [Post] = PlotPosterior(filename, Prior, Theory, d)
    [depth,gamma,phi,rho,caliper] = ImportWellLog(filename);

    x1 = 3.0:0.05:8.0;
    x2 = 2.0:0.02:3.0;
    x3 = 0.0:0.01:0.5;
    [X, Y, Z, Post, PostP] = CalculatePosterior(Prior, Theory, x1(:), x2(:), x3(:));
    
    figure;
    set(slice(X,Y,Z,Post,3.0:.5:8.0,[],[]),'edgecolor','none');
    colorbar;
    title(['Posterior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);

    figure;
    set(slice(X,Y,Z,PostP,3.0,3.0,0.0),'edgecolor','none');
    colorbar;
    title(['Posterior Joint PDF at depth of ',num2str(depth(d)),' ft']);
    xlabel('Velocity (km/s)');
    ylabel('Density (g/cc)');
    zlabel('Porosity (.)');
    axis tight;
    daspect([1 1 0.33]);
    view(32,32);
    
    [m,Prx] = max(Prior(:));
    [I1,I2,I3] = ind2sub(size(Prior),Prx);
    
    [m,Psx] = max(Post(:));
    [II1,II2,II3] = ind2sub(size(Post),Psx);
    
    PrV = Prior(I1,:,I3);
    PrR = reshape(Prior(:,I2,I3), 1, length(x2));
    PrP = reshape(Prior(I1,I2,:), 1, length(x3));
    stdPr = (std(x1,PrV).^2+std(x2,PrR).^2+std(x3,PrP).^2).^(0.5);
    
    PsV = Post(II1,:,II3);
    PsR = reshape(Post(:,II2,II3), 1, length(x2));
    PsP = reshape(Post(II1,II2,:), 1, length(x3));
    stdPs = (std(x1,PsV).^2+std(x2,PsR).^2+std(x3,PsP).^2).^(0.5);
    PRPDF = normpdf(x1,x1(I2),stdPr);
    PSTPDF = normpdf(x1,x1(II2),stdPs);
    
    figure;
    plot(x1,PRPDF, 'b'); hold on;
    plot(x1, PSTPDF, 'r'); hold on;
    plot([x1(I2)-std(x1,Prior(I1,:,I3)) x1(I2)+std(x1,Prior(I1,:,I3))], [max(PRPDF) max(PRPDF)], 'b'); hold on;
    plot([x1(I2) x1(I2)], [0 max(PRPDF)], 'b'); hold on;
    plot([x1(II2)-std(x1,Post(II1,:,II3)) x1(II2)+std(x1,Post(II1,:,II3))], [max(PSTPDF) max(PSTPDF)], 'r'); hold on;
    plot([x1(II2) x1(II2)], [0 max(PSTPDF)], 'r'); hold on;
    xlabel('v (km/s)');
    legend('Prior','Posterior');
    title(['Prior vs. Posterior at ',num2str(depth(d)),' ft']);
end

function LogPanel(filename)
    [depth,gamma,phi,rho,caliper] = ImportWellLog(filename);
    
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

