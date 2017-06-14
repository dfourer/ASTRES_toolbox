function Imin=find_minima(x)

[Nx,Ny]=size(x);

Imin=zeros(Nx,Ny);

%lines
for j=1:Nx
  %Imin(j,:) = Imin(j,:) + imregionalmin(x(j,:));
  Imin(j,:) = Imin(j,:) + find_minima_line(x(j,:));
end

%columns
for j=1:Ny
  %Imin(:,j) = Imin(:,j) + imregionalmin(x(:,j)); 
  Imin(:,j) = Imin(:,j) + find_minima_column(x(:,j)); 
end

%first diagonal
for j=-(Nx-1):(Ny-1)
   diagonal = diag(x,j);
   %diagmin = imregionalmin(diagonal);
   diagonal=transpose(diagonal);
   diagmin = find_minima_line(diagonal);
   if (Nx == Ny)
       Imin = Imin + diag(diagmin,j);
   else Imdiag = diag(diagmin,j); %case of non-square images
       [Nxdiag,Nydiag] = size(Imdiag);
       for abs=1:min(Nx,Nxdiag)
           for ord=1:min(Ny,Nydiag)
               Imin(abs,ord)=Imin(abs,ord)+Imdiag(abs,ord);
           end
       end
   end
end

%second diagonal

xflip = flipdim(x,2);
for j=-(Nx-1):(Ny-1)
   diagonal = diag(xflip,j);
   %diagmin = imregionalmin(diagonal);
   diagonal=transpose(diagonal);
   diagmin = find_minima_line(diagonal);
   Iminflip = diag(diagmin,j);
   if (Nx == Ny)
       Imin = Imin + flipdim(Iminflip,2);
   else Imflip = flipdim(Imdiag,2);
       [Nxdiag,Nydiag] = size(Imdiag);
       for abs=1:min(Nx,Nxdiag) %case of non-square images
           for ord=1:min(Ny,Nydiag)
               Imin(abs,ord)=Imin(abs,ord)+Imflip(abs,ord);
           end
       end
   end
   %Imin = Imin + flipdim(Iminflip,2);
end

ind = find (Imin > 1);
Imin(ind)=1;

end