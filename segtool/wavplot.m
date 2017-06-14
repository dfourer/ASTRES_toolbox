function wavplot(wtype,p1,p2)
% wavplot : display the wavelet wtype in both time and frequency domains

switch(wtype)
    case 'cmor' % Complex Morlet
        N = 16*1024;
        t = linspace(-5*sqrt(p1),5*sqrt(p1),N);
        nu = linspace(p2-5/pi/sqrt(p1),p2+5/pi/sqrt(p1),N);
        figure();
        subplot(1,2,1);
        psi = 1/sqrt(p1*pi)*exp(2*1i*pi*p2*t) .* exp(-t.^2/p1);
        plot(t,real(psi),'b--',t,imag(psi),'r-.',t,abs(psi),'k-');
        legend('real part','imaginary part','modulus');
        xlabel('t');ylabel('\psi(t)');
        subplot(1,2,2);
        psi = cmor(p1,p2,nu);
        plot(nu,psi);
        xlabel('\xi');ylabel('\Psi (\xi)');
        %set(gca,'xlim',[p2-0.2*p1,p2+0.2*p1]);
    case 'gmor' % Generalized Morse
        N = 16*1024;
        oml = (p1/p2)^(1/p2);
        nu = linspace(-1000/sqrt(p1*p2),1000/sqrt(p1*p2),N);
        t = (-N/2:N/2-1)/(nu(end)-nu(1));
        figure();
        subplot(1,2,2);
        psi = gmor(p1,p2,nu);
        plot(nu,psi);
        xlabel('\xi');ylabel('\Psi(\xi)');
        set(gca,'xlim',[oml-10/sqrt(p1*p2),oml+10/sqrt(p1*p2)]);
        
        subplot(1,2,1);
        psir = sqrt(N)*ifftshift(ifft(psi)).*(-1).^(0:N-1);
        plot(t,real(psir),'b--',t,imag(psir),'r-.',t,abs(psir),'k-');
        legend('real part','imaginary part','modulus');
        xlabel('t');ylabel('\psi(t)');
        set(gca,'xlim',[-10*sqrt(p1*p2),10*sqrt(p1*p2)]);
    case 'bump' % Bump Wavelet
        N = 16*1024;
        nu = linspace(p1-2*p2/p2,p1+2*p2/p2,N);
        t = (-N/2:N/2-1)/(nu(end)-nu(1));
        figure();
        subplot(1,2,2);
        psi = bump(p1,p2,nu);
        plot(nu,psi);
        xlabel('\xi');ylabel('\Psi(\xi)');
        set(gca,'xlim',[p1-2*p2,p1+2*p2]);
        
        subplot(1,2,1);
        psir = sqrt(N)*ifftshift(ifft(psi)).*(-1).^(0:N-1);
        plot(t,real(psir),'b--',t,imag(psir),'r-.',t,abs(psir),'k-');
        legend('real part','imaginary part','modulus');
        xlabel('t');ylabel('\psi(t)');
        set(gca,'xlim',5./[-p2,p2]);
        
end

end
        