function [imamplitude,imorientation,imfreq,immodel]=sinus_estimation_parallel(t1,N,sizex,sizey)

indx=1;
indy=1;
for indexi = floor(N/2):N:sizex-floor(N/2)
    for indexj = floor(N/2):N:sizey-floor(N/2)
        c=recover_sinus(t1(:,:,indx,indy));
        
        %Compute instantaneous frequencies
        %Reshape t1 and t2 into 2N2(N1-2)X3 matrix
        ct=transpose(c);
        t1r1=zeros(N*(N-2),3);
        for indi=1:N
            for indk=1:N-2
                t1r1(indk+(N-2)*(indi-1),:) = ct(indi,N-indk-1:N-indk+1);
            end
        end
        t1r2=rot90(t1r1,2);
        t1r=[t1r1 ; t1r2];

        t2r1=zeros(N*(N-2),3);
        for indi=1:N
            for indk=1:N-2
                t2r1(indk+(N-2)*(indi-1),:) = c(indi,N-indk-1:N-indk+1);
            end
        end
        t2r2=rot90(t2r1,2);
        t2r=[t2r1 ; t2r2];

        %Compute annihilating filters
        [U1,S1,V1]=svd(t1r);
        [U2,S2,V2]=svd(t2r);
        h1=V1(:,3);
        h2=V2(:,3);

        %roots of annihilating filters
        roots1=roots(h1);
        roots2=roots(h2);

        %frequency estimation
        freq1=angle(roots1(1))/(2*pi);
        freq2=angle(roots2(1))/(2*pi);
        freq=sqrt(freq1^2 + freq2^2);

        %disambigaute frequency and compute amplitude and orientation
        cana=hilbert(c);

        cplus=zeros(N,N);
        cmoins=zeros(N,N);

        for indi=1:N
            for indj=1:N
                cplus(indi,indj)=exp(1i*2*pi*(freq1*indi+freq2*indj));
                cmoins(indi,indj)=exp(1i*2*pi*(freq1*indi-freq2*indj));
            end
        end
        canaplus=cana./cplus;
        canamoins=cana./cmoins;

        if var(canaplus(:)) <= var(canamoins(:))
            amplitude = mean(abs(canaplus(:)));
            if freq1 == 0 
                orientation = pi/2;
            else orientation = atan(freq2/freq1);
            end
        else amplitude = mean(abs(canamoins(:)));
            if freq1 == 0
                orientation = -pi/2;
            else orientation = atan(-freq2/freq1);
            end
        end
        
        immodel(indexi-floor(N/2)+1:indexi+floor(N/2),indexj-floor(N/2)+1:indexj+floor(N/2))=c(:,:);
        imamplitude(indexi-floor(N/2)+1:indexi+floor(N/2),indexj-floor(N/2)+1:indexj+floor(N/2))=amplitude;
        imfreq(indexi-floor(N/2)+1:indexi+floor(N/2),indexj-floor(N/2)+1:indexj+floor(N/2))=freq;
        imorientation(indexi-floor(N/2)+1:indexi+floor(N/2),indexj-floor(N/2)+1:indexj+floor(N/2))=orientation;
        indy=indy+1;
    end
    indx=indx+1;
    indy=1;
end

end