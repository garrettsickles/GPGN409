function Laplacian2DFilterTest()
    m_dim = 13;
    d_dim = 15;
    k_dim = m_dim*d_dim;
    
    tm = zeros(m_dim,d_dim);
    tm(1,:) = [1 0 0 0 0 0 0 0 0 1 0 0 0 0 1];
    tm(5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    tm(7,:) = [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
    tm(9,:) = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    tm(13,:) = [1 0 0 0 0 1 0 0 0 0 0 0 0 0 1];
    
    B_lap = LaplacianFilter2D(m_dim, d_dim);
    
    B_lap2 = zeros(k_dim);
    for i=1:k_dim
        for j=1:k_dim
            B_lap2(i,j) = LaplacianFilter2DAt(i,j,m_dim,d_dim);
        end
    end
            
    figure;
    imshow(mat2gray(tm),...
        'InitialMagnification', 3000, 'Colormap', flipud(gray));
    set(gca,'Ydir','reverse');

    figure;
    imshow(mat2gray(reshape(B_lap*reshape(tm,k_dim,1),m_dim,d_dim)),...
        'InitialMagnification', 3000);
    colormap(gray)
    set(gca,'Ydir','reverse');
    
    figure;
    imshow(mat2gray(reshape(B_lap2*reshape(tm,k_dim,1),m_dim,d_dim)),...
        'InitialMagnification', 3000);
    colormap(gray)
    set(gca,'Ydir','reverse');
    
    figure;
    imshow(mat2gray(B_lap),...
        'InitialMagnification', 300);
    
    figure;
    imshow(mat2gray(B_lap2),...
        'InitialMagnification', 300);
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