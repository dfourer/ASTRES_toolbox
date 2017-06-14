function display_frequency(frequency,coherency)
    %display frequencies
    figure;
    imagesc(frequency,[0 1]);
    colormap('hsv');
    axis off;
    axis image;
    colorbar;
    title(sprintf('Frequency of IMF'));
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
    hsv = zeros(size(frequency, 1), size(frequency, 2), 3);
    hsv(:,:,1) = frequency;
    hsv(:,:,2) = coherency;
    hsv(:,:,3) = ones(size(hsv,1), size(hsv,2));
    rgb = hsv2rgb(hsv);
    imagesc(rgb)
    axis off;
    axis image;
    title(sprintf('Composite frequency and coherency of IMF'));
    clear hsv