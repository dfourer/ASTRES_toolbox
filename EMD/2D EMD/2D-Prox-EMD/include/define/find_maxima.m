function Imax=find_maxima(x)

[Nx,Ny]=size(x);

Imax=zeros(Nx,Ny);

%lines
for j=1:Nx
  %Imax(j,:) = Imax(j,:) + imregionalmax(x(j,:));  
  Imax(j,:) = Imax(j,:) + find_maxima_line(x(j,:));  
end

%columns
for j=1:Ny
  %Imax(:,j) = Imax(:,j) + imregionalmax(x(:,j));
  Imax(:,j) = Imax(:,j) + find_maxima_column(x(:,j)); 
end

%first diagonal
for j=-(Nx-1):(Ny-1)
   diagonal = diag(x,j);
   %diagmax = imregionalmax(diagonal);
   diagonal=transpose(diagonal);
   diagmax=find_maxima_line(diagonal);
   if (Nx == Ny)
       Imax = Imax + diag(diagmax,j);
   else Imdiag = diag(diagmax,j); %case of non-square images
       [Nxdiag,Nydiag] = size(Imdiag);
       for abs=1:min(Nx,Nxdiag)
           for ord=1:min(Ny,Nydiag)
               Imax(abs,ord)=Imax(abs,ord)+Imdiag(abs,ord);
           end
       end
   end
end

%second diagonal

xflip = flipdim(x,2);
for j=-(Nx-1):(Ny-1)
   diagonal = diag(xflip,j);
   %diagmax = imregionalmax(diagonal);
   diagonal=transpose(diagonal);
   diagmax=find_maxima_line(diagonal);
   Imaxflip = diag(diagmax,j);
   if (Nx == Ny)
       Imax = Imax + flipdim(Imaxflip,2);
   else Imflip = flipdim(Imdiag,2);
       [Nxdiag,Nydiag] = size(Imdiag);
       for abs=1:min(Nx,Nxdiag) %case of non-square images
           for ord=1:min(Ny,Nydiag)
               Imax(abs,ord)=Imax(abs,ord)+Imflip(abs,ord);
           end
       end
   end
   %Imax = Imax + flipdim(Imaxflip,2);
end

ind = find (Imax > 1);
Imax(ind)=1;

end


% function linemax=find_maxima_line(line)
% 
% [~,sizel]=size(line);
% linemax=zeros(1,sizel);
% for l=2:sizel-1
%     if (line(l)>line(l-1) & line(l)>line(l+1))
%         linemax(l)=1;
%     end
% end
% end
% 
% function columnmax=find_maxima_column(column)
% 
% [sizec,~]=size(column);
% columnmax=zeros(sizec,1);
% for c=2:sizec-1
%     if (column(c)>column(c-1) & column(c)>column(c+1))
%         columnmax(c)=1;
%     end
% end
% end

% function fill_diagonal (x,diag,j)
% 
% [Nx,Ny]=size(x);
% 
% if (j == 0)
%     for k=1:min(Nx,Ny)
%         x(k,k) = diag(k);
%     end
% elseif (j<0)
%     
% elseif (j>0)
%     
% end