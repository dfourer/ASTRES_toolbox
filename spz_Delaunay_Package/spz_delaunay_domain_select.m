function [Domains, Domains_Tri] = spz_delaunay_domain_select(triZ, xZ, yZ, Nx, ec1, ec2, Pruning_Thresh)
%  [Domains, Domains_Tri] = spz_delaunay_domain_select(triZ, xZ, yZ, Nx, ec1, ec2, Pruning_Thresh)
%   Extract domains (Domains) as a set of triangles (Domains_Tri) via Delaunay's approach
%
% P. Flandrin & Ph. Depalle
% 2015, July 1st
%
% inputs
%   triZ           = Vector of triangles
%   xZ             = Vector of X vertices' coordinates
%   yZ             = Vector of Y vertices' coordinates
%   Nx             = Block Size (has to be a power of 2, and preferably <= 1024)
%   ec1            = threshold for large edges' length (ref = 2)
%   ec2            = threshold for small edges' length (ref = 0.25)
%   Pruning_Thresh = minimum size for a given domain (e.g. 2)
%                  
% outputs
%   Domains        = Return array cell of selected domains.
%   Domains_Tri    = Return array of triangle involved in selected domains (including the prunned ones).
% Calls:
%   spz_delaunay_dom4b2
% 	spz_delaunay_domain_pruning


fprintf('\nSelection.');
[Domains_All, N_Domains_All, LM, Domains_Tri]...
                        = spz_delaunay_dom4b2(Nx, triZ, xZ, yZ, ec1, ec2);

%     Root_File_Name_Sel = [Root_File_Name '_Sel' '_ec1_' num2str(ec1) '_ec2_' num2str(ec2)];
%     spz_delaunay_save(triZLM, xZ, yZ, Domains_All, Nx, Root_File_Name_Sel);

if(Pruning_Thresh == 0)
    Domains     = Domains_All;
    fprintf('\n\tNo Pruning\n');
else
    Domains     = spz_delaunay_domain_pruning(Domains_All, Pruning_Thresh);
    fprintf('\n\tPruning\n');
end
end