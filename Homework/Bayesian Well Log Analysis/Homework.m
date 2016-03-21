function Homework(filename)
    [depth,gamma,phi,rho,caliper] = Import_Well_Log(filename);
    
    figure;
    v = Vel_From_Gamma(gamma);
    x1 = 0.0:0.01:4.0;
    x2 = 0.0:0.01:0.3;
    x3 = 3.0:0.01:7.5;
    [X,Y,Z] = meshgrid(x1,x2,x3);
    sigma = diag(ones(3,1)*abs(caliper(351)-6.0));
    % [0.05 0 0;0 0.001 0; 0 0 abs(caliper(351)-6.0)]
    F = mvnpdf([X(:) Y(:) Z(:)], [rho(351) phi(351) v(351)], sigma);
    F = reshape(F,length(x2),length(x1),length(x3));
    xslice = [1.2,1.8,3.0]; 
    yslice = [0.15, 0.20]; 
    zslice = [3.5, 5.5];
	set(slice(X,Y,Z,F,xslice,yslice,zslice),'edgecolor','none');
    
    figure;
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

function [vp] = Vel_From_Wyllie(phi,vm,vf)
    vp = ((1-phi)./(vm)+phi./vf).^(-1);
end

function [rho] = Rho_From_Phi(phi,rm,rf)
    rho = ((1-phi).*(rm)+phi.*rf);
end

function [vp] = Vel_From_Gamma(gamma)
    vp = 5.654 - 0.008.*gamma;
end

%function Prior()