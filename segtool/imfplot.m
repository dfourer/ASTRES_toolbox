function imfplot(imf,dt,ha)
% imfplot displays the IMFs on a same graphic

mstyl = {'b','r','k','g','y','b','r','k','g','y','b','r','k','g','y','b','r','k','g','y'};
t = dt*(0:size(imf,2)-1);
leg = {};

cla(ha);hold(ha,'on');

for j=1:size(imf,1)
    plot(t,imf(j,:),mstyl{j},'Parent',ha);
    leg = [leg {['IMF ' num2str(j)]}];
end

legend(ha,leg);