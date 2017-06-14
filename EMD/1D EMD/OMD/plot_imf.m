function plot_imf(imf,t,suffixe)

if nargin<2
    t = linspace(0,1,size(imf,2));
    suffixe = '';
elseif nargin<3
    suffixe = '';
end

M = size(imf,1);

if abs(t-1)<=eps
    t = 1:size(imf,2);
end

subplot(M,1,1);
hold on;
title(['The decomposition by ' suffixe]);

for i=1:M
    if i>1
        subplot(M,1,i);
    end
    plot(t,imf(i,:));
    xlabel('time');
    ylabel(['IMF ' num2str(i)]);
end

end