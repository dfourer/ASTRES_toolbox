function p = prox_L12(x,gamma)


[n,m,~] = size(x);
wx = zeros(2,n*m);
wx(1,:) = reshape(x(:,:,1),1,n*m);
wx(2,:) = reshape(x(:,:,2),1,n*m);
wp = zeros(size(wx));
nb_r = size(wp,1);
tmp = sqrt(sum(wx.^2,1));
ind = find(tmp>gamma);
wp(:,ind) = kron(ones(nb_r,1),(1 - gamma./tmp(ind))).*wx(:,ind);
p = zeros(size(x));
p(:,:,1) = reshape(wp(1,:),n,m);
p(:,:,2) = reshape(wp(2,:),n,m);

