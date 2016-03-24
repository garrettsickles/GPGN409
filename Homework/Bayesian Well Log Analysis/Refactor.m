function Refactor()
    filename = 'data.txt';
    d = 33;
    
    Prior = PlotPrior(filename, d);
    Theory = PlotTheory();
    PlotPosterior(filename, Prior, Theory, d);
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
    sigma = [std(vel(:)) 0 0;0 std(rho(:)) 0; 0 0 std(phi(:))].*((abs(caliper(d)-5.999))/range(caliper(:)));
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
    
    PRPDF = normpdf(x1,x1(I2),std(x1,Prior(I1,:,I3)));
    PSTPDF = normpdf(x1,x1(II2),std(x1,Post(II1,:,II3)));
    
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