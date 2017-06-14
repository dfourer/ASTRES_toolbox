function wp = prox_L1(wx,gamma)

wp = max(abs(wx)-gamma,0).*sign(wx);


