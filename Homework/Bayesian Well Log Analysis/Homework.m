function Homework()
    filename = 'data.txt';
    PlotPrior(filename, 350);
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

function PlotPrior(filename, d)
    [depth,gamma,phi,rho,caliper] = Import_Well_Log(filename);
    
    figure;
    v = Vel_From_Gamma(gamma);
    x1 = 3.0:0.005:7.5;
    x2 = 2.0:0.01:3.0;
    x3 = 0.0:0.001:0.3;
    [X,Y,Z] = meshgrid(x1,x2,x3);
    sigma = [range(v(:))*2 0 0;0 range(rho(:)) 0; 0 0 0.1].*((abs(caliper(d)-6.0))/range(caliper(:)));
    range(caliper(:))
    F = mvnpdf([X(:) Y(:) Z(:)], [v(d) rho(d) phi(d)], sigma);
    F = reshape(F,length(x2),length(x1),length(x3));
	set(slice(X,Y,Z,F,[3.0:0.5:7.5],[],[]),'edgecolor','none');
    axis tight;
    daspect([1 1 0.33]);
    view(30,32);
    figure;
    FF = F;
    F(length(x2),:,:) = sum(FF(:,:,:),1);
    F(:,1,:) = sum(FF(:,:,:),2);
    F(:,:,1) = sum(FF(:,:,:),3);
	set(slice(X,Y,Z,F,[3.0],[3.0],[0.0]),'edgecolor','none');
    axis tight;
    daspect([1 1 0.33]);
    view(30,32);
end