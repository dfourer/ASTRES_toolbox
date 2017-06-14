function display_orientation(ang,coherency)
    %display angles
    figure;
    imagesc(ang,[-pi/2 pi/2]);
    colormap('hsv');
    colorbar;
    axis off;
    axis image;
    title(sprintf('Angle of IMF'));
    %display coherency
    figure;
    imagesc(coherency,[0 1]);
    colormap('gray');
    colorbar;
    axis off;
    axis image;
    title(sprintf('Coherency of IMF'));
    %display a composite image
    figure;
    clear hsv;
    hsv = zeros(size(ang, 1), size(ang, 2), 3);
    hsv(:,:,1) = (ang+pi/2)/pi;
    hsv(:,:,2) = coherency;
    hsv(:,:,3) = ones(size(hsv,1), size(hsv,2));
    rgb = hsv2rgb(hsv);
    imagesc(rgb)
    axis off;
    axis image;
    title(sprintf('Composite angle and coherency of IMF'));
    clear hsv