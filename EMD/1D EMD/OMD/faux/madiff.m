function deriv = madiff(s,t,k)
%% deriv : gives the k-order derivative. Completion with constants. Uses a 4-th order formula.

mykern = [-1 8 0 -8 1];
deriv = s;

for i=1:k
    aux = conv(deriv,mykern);
    deriv = aux;
    deriv = deriv/(12*(t(2)-t(1)));
end

deriv=deriv(2*k+1:end-2*k);
deriv(1:2*k) = deriv(2*k+1);
deriv(end-2*k+1:end) = deriv(2*k);


end