function wp = proj_L1_modifEMD(wx,gamma)

[n,m] = size(wx);
wx = wx(:)';
wp = NormL1_project(wx,1,gamma);
wp = reshape(wp,n,m);