function [List_of_Segments, Nsegments, LM, triZLM] = spz_delaunay_dom4b2(Nx, triZ, xZ, yZ, ec, ec2)
%
% P. Flandrin & Ph. Depalle, 
% 2015, June, 10th
%
% 
% input
%   Signal_Length  Length of Signal
%   Nx             Block Size (has to be a power of 2, and preferably <= 1024)
%   triZ           Set of triangles
%   xZ, yZ         Coordinates of triangles vertices
%   ec             threshold for large edges' length (ref = 2)
%   ec2            threshold for small edges' length (ref = 0.25)
%   fignum         figure number for colored sub-domains
%
% output
%   inseg          masks
%   Nsegments      Number of sub-domains
%   LM             Vector of Maximum Length per Triangle
%
% calls
%   roundgauss.m

global gr2
global C
%%  
  
Nfft = 2*Nx ;%2*Nx ; % number of FFT bins
prec = 10^(-6) ;
[H L] = roundgauss(Nfft, prec);
Ly = Nfft/L ;
Lx = L ;
LM = [] ; % max distances between min

 for kk = 1:length(triZ)
% %    kk = jnd(k) ;
     L1 = sqrt((xZ(triZ(kk,1))-xZ(triZ(kk,2)))^2/Lx^2+(yZ(triZ(kk,1))-yZ(triZ(kk,2)))^2/Ly^2) ;
     L2 = sqrt((xZ(triZ(kk,1))-xZ(triZ(kk,3)))^2/Lx^2+(yZ(triZ(kk,1))-yZ(triZ(kk,3)))^2/Ly^2) ;
     L3 = sqrt((xZ(triZ(kk,3))-xZ(triZ(kk,2)))^2/Lx^2+(yZ(triZ(kk,3))-yZ(triZ(kk,2)))^2/Ly^2) ;
     LM = [LM max(max(L1,L2),L3)] ;
end

% threshold
indLZM = find(LM>ec) ;
%triZLM = triZ(jnd(indLZM),:) ;
triZLM = triZ(indLZM,:) ;
LM = LM(indLZM);

fprintf('\nFind common segments.');
% mask from adjacent triangles
[AT,BT] = size(triZLM) ;
conc = [] ;
gr = [] ;
p = 1 ;
for k = 1:AT-1
    for m = k+1:AT
 %        Sis = sum(ismember(triZLM(k,:),triZLM(m,:))) ; 
         [logic_vertex, add_vertex] = ismember(triZLM(k,:),triZLM(m,:)) ;
        if (sum(logic_vertex)==2)
            add_vertex = find(add_vertex); % Remove '0' from the list
            common_length =...
                sqrt( (xZ(triZLM(m,add_vertex(1))) - xZ(triZLM(m,add_vertex(2))))^2/Lx^2 +...
                      (yZ(triZLM(m,add_vertex(1))) - yZ(triZLM(m,add_vertex(2))))^2/Ly^2);
            if(common_length > ec2) %ec/128)
                gr(p,:) = [k m];
                conc = [conc' [k m]']' ; % Might be useless
                p = p+1 ;
            end
        end
    end
end
gr2 = gr;
fprintf('\nFind domains.');
U = unique(conc)' ;
LU = length(U) ;
[Ag,Bg] = size(gr) ;
C = [] ;
Activation_Flag = false;
%Find_Domains(gr);
for k = 1:Ag-1
    fprintf('\n\tFind domains: Segment %d over %d', k, Ag);
    grtest = gr(k,:) ;
    grtest2 = grtest ;
    dgr = 1 ;
    while dgr~=0;
        for m = k+1:Ag
            ugr = union(grtest2,gr(m,:)) ;
            if length(ugr)<length(grtest)+length(gr(m,:)) & gr(m,:)~=[0 0] ; 
                grtest2 = ugr ;
                gr(m,:) = [0 0] ;
                Activation_Flag = true;
            end
        end
        dgr = length(grtest) - length(grtest2) ;
        grtest = grtest2 ;
    end
    if(Activation_Flag)
        C = [C grtest];  % Vector of domain
    else
        C = [C grtest 0];
        grtest
    end
end

fprintf('\n\nBuild domains.');

indC = find(C==0) ;
Nsegments = min(find(diff(indC)==1)) ;
%intot = zeros(NAx,NAy) ;

for seg = 1:Nsegments
    if seg==1 ;
        Useg = C(1:indC(seg)-1) ;
    else
        Useg = C(indC(seg-1)+1:indC(seg)-1);   
    end
    fprintf('\n\tBuild Domain: %d out of %d. Last segment: %d', seg, Nsegments, Useg(end));
    List_of_Segments{seg} = [Useg]; 
end

end




