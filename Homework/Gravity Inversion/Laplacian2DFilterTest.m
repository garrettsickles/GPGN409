function Laplacian2DFilterTest()
    m_dim = 13;
    d_dim = 15;
    k_dim = m_dim*d_dim;
    
    control = zeros(m_dim,d_dim);
    control(1,:) = [-4 1 0 0 0  0 0  0 1 -4 1 0 0 1 -4];
    control(2,:) = [ 1 0 0 0 0  0 0  0 0  1 0 0 0 0  1];
    control(3,:) = [ 0 0 0 0 0  0 0  0 0  0 0 0 0 0  0];
    control(4,:) = [ 0 0 0 0 0  0 0  0 0  0 0 0 0 0  1];
    control(5,:) = [ 0 0 0 0 0  0 0  0 0  0 0 0 0 1 -4];
    control(6,:) = [ 0 0 0 0 0  0 0  1 0  0 0 0 0 0  1];
    control(7,:) = [ 0 0 0 0 0  0 1 -4 1  0 0 0 0 0  0];
    control(8,:) = [ 1 0 0 0 0  0 0  1 0  0 0 0 0 0  0];
    control(9,:) = [-4 1 0 0 0  0 0  0 0  0 0 0 0 0  0];
    control(10,:) = [ 1 0 0 0 0  0 0  0 0  0 0 0 0 0  0];
    control(11,:) = [ 0 0 0 0 0  0 0  0 0  0 0 0 0 0  0];
    control(12,:) = [ 1 0 0 0 0  1 0  0 0  0 0 0 0 0  1];
    control(13,:) = [-4 1 0 0 1 -4 1  0 0  0 0 0 0 1 -4];

    
    % tm = randi([0 1],m_dim,d_dim);
    tm = zeros(m_dim,d_dim);
    tm(1,:)  = [1 0 0 0 0 0 0 0 0 1 0 0 0 0 1];
    tm(5,:)  = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    tm(7,:)  = [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
    tm(9,:)  = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    tm(13,:) = [1 0 0 0 0 1 0 0 0 0 0 0 0 0 1];
    
    B_lap = LaplacianFilter2D(m_dim, d_dim);
    B_lap2 = zeros(k_dim);
    
    for i=1:k_dim
        for j=1:k_dim
            B_lap2(i,j) = LaplacianFilter2DAt(i,j,m_dim,d_dim);
        end
    end
    
    m1 = reshape(B_lap*reshape(tm,k_dim,1),m_dim,d_dim);
    m2 = reshape(B_lap2*reshape(tm,k_dim,1),m_dim,d_dim);
    figure;
    imshow(mat2gray(tm),...
        'InitialMagnification', 3000, 'Colormap', flipud(gray));
    title('Starting Model');
    set(gca,'Ydir','reverse');

    figure;
    subplot(2,3,1);
    imshow(mat2gray(control),...
        'InitialMagnification', 3000);
    title('Control');
    colormap(gray);
    set(gca,'Ydir','reverse');
    
    subplot(2,3,2);
    imshow(mat2gray(m1),...
        'InitialMagnification', 3000);
    title('Matrix');
    colormap(gray);
    set(gca,'Ydir','reverse');
    
    subplot(2,3,3);
    imshow(mat2gray(m2),...
        'InitialMagnification', 3000);
    title('Function');
    colormap(gray);
    set(gca,'Ydir','reverse');
    
    subplot(2,3,5);
    imshow(mat2gray(control-m1),...
        'InitialMagnification', 3000);
    title('Difference');
    colormap(flipud(gray));
    set(gca,'Ydir','reverse');
    
    subplot(2,3,6);
    imshow(mat2gray(control-m2),...
        'InitialMagnification', 3000);
    title('Difference');
    colormap(flipud(gray));
    set(gca,'Ydir','reverse');
    
    figure;
    subplot(2,2,1);
    imshow(mat2gray(B_lap),...
        'InitialMagnification', 300);
    title('Matrix');
    set(gca,'Ydir','reverse');
    
    subplot(2,2,2);
    imshow(mat2gray(B_lap2),...
        'InitialMagnification', 300);
    title('Function');
    set(gca,'Ydir','reverse');
    
    subplot(2,2,[3 4]);
    imshow(mat2gray(B_lap2-B_lap),...
        'InitialMagnification', 300);
    title('Difference');
    colormap(flipud(gray));
    set(gca,'Ydir','reverse');
end

function [ B ] = LaplacianFilter2D(x,y)
    k = x*y;
    B = zeros(k);
    B = B + diag(ones(k,1).*(-4));
    B = B + diag(ones(k-1,1),1);
    B = B + diag(ones(k-1,1),-1);
    B = B + diag(ones(k-x,1),x);
    B = B + diag(ones(k-x,1),-x);
    for i=1:(y-1)
        j = i*x; 
        B(j,j+1) = 0;
        B(j+1,j) = 0;
    end
end

function [ v ] = LaplacianFilter2DAt(r, c, x, y)
    if(r == c)
        v = -4;
    elseif(r-1 == c) % Below or to the left
        if(mod(r-1,x) == 0)
            v = 0;
        else
            v = 1;
        end
    elseif(c-1 == r) % Above or to the right
        if(mod(c-1,x) == 0)
            v = 0;
        else
            v = 1;
        end
    elseif(c-x == r || c+x == r)
        v = 1;
    else
        v = 0;
    end
end