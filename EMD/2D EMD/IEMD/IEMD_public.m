function  [ix,resx,medel]=IEMD_public(image,epsilon,numberimfs,conn)
% Decompose an image into its intrinsic modes IMFs and corresponding
% residues.
%  
%
% INPUT:
% image- a matlab 2D vector;
% epsilon- value of stop criterion, between 0 and 1 but you can try other
% and see what happens
% numberimfs- maximum number of IMF
% conn- Type of neighbourhood in extrema point selection
%  '4m'     two-dimensional four-connected neighborhood, use 'imregionalmax' in Images toolbox
%  '8m'     two-dimensional eight-connected neighborhood, use 'imregionalmax' in Images toolbox
%  '4'     two-dimensional four-connected neighborhood, simple implementation
%  '8'     two-dimensional eight-connected neighborhood, simple implementation
% Note all that these give different results 
%
% OUTPUT:
% ix - 3D vector with all IMFs
% resx - 3D vector with all residues
%
% 
%
% AUTHOR: 
% Anna Linderhed 
% 
% DATE:
% 2001-07-20 -- created
% last update 2005-10-10 
 
%
% NOTE:
% The interpolation use 'tpaps' in the Matlab 6.5 Splines toolbox. This is
% extremely slow and memory consuming; unless you have a very good computer
% do not use this code on larger images than 128x128. Please make a
% better implementation and tell me about it.

% This code may be used for education and inspiration. If using this code results in published work
%  a reference to the source must be included in the paper:
% http://go.to/imageemd 

%%


%%
tic
ix=zeros(size(image,1),size(image,2),numberimfs);
resx=zeros(size(image,1),size(image,2),numberimfs);
f=0;

[imf,res,medel]=findimf(image,1,epsilon,conn);

for i=1:numberimfs
   resx(:,:,i)=res;
   ix(:,:,i)=imf;
    
   if i==numberimfs
      f=1;
   end
   [max,min]=maxomin(res,conn);
   [a dummy dummy]=find(max);
   if size(a)<1
      f=1;
   end
   
   [a2 dummy dummy]=find(min);
   if size(a2)<1
      f=1;
   end


   if f==1;
      i  
      break
   end
   
  [imf,res,medel]=findimf(res,1,epsilon,conn);
end 
toc
  




%%*******************************************************************
%%

function [imf,res,medel]=findimf(xin,p,epsilon,conn)

sis=size(xin,1);
numpix=sis*sis;
hk=xin-1;
medel=hk;
k=0;
xin1=xin;
for i=1:100
    
  if      sum(sum(abs(medel)))/numpix> epsilon
      
      % k=sum(sum(abs(medel)))/numpix;
       xin=hk; 
       [hk, medel,f]=interpolateimf(xin,p,conn);
      % b=i;

  end

imf=hk;
res=xin1-hk;

end
%figure;imagesc(imf);colormap(gray);axis off;

%%%%********************************************************************
%%
function [c,medel,f]=interpolateimf(x,p,conn)

[y,v]=maxomin(x,conn);
      
 [Iy,Jy,Vy] = find(y);
[Iv,Jv,Vv] = find(v);
[Iy,Jy,Vy] = find(y);
  xyy=[Iy';Jy'];
  xyv=[Iv';Jv'];



[m1,n1] = size(y);[m2,n2] = size(v);
[m12,n12,v12] = find(Iy);
[m22,n22,v22] = find(Iv);

   
if size(Iy)<200
  for k=1:2:size(y,1)
    y(1,k)=x(1,k);
    y(size(x,1),k)=x(size(x,1),k);
    y(k,size(x,2))=x(k,size(x,2));
       y(k,1)=x(k,1);
  end
end
[Iv,Jv,Vv] = find(v);
if size(Iv)<200
  for k=1:2:size(v,1)
    v(1,k)=x(1,k);
    v(size(x,1),k)=x(size(x,1),k);
    v(k,size(x,2))=x(k,size(x,2));
       v(k,1)=x(k,1);
  end
end
        
        


if max(size(v12,1),size(v12,2))<3
     medel=zeros(size(x));
     f=1;
c=x-medel;
    return
end
if max(size(v22,1),size(v22,2))<3
   medel=zeros(size(x));
c=x-medel;
f=1;
    return
end



sty = tpaps(xyy,Vy',p); %figure(1);fnplt(sty)
 valy=fnval(sty,{1:m1,1:n1});
 stv = tpaps(xyv,Vv',p); %figure(2);fnplt(stv)
 valv=fnval(stv,{1:m2,1:n2});
f=0;
medel=(valy+valv)/2;
%figure(1);imagesc(valy);colormap(gray);axis off;
%figure(2);imagesc(valv);colormap(gray);axis off;
%figure(3);imagesc(medel);colormap(gray);axis off;
c=x-medel;
%figure(4);imagesc(c);colormap(gray);axis off;





%**************************
%%
%*************************************************************************


function [max,min]=maxomin(x,conn)

y=zeros(size(x));

v=zeros(size(x));



switch conn
    case '4'
 
  
for j=2:size(x,2)-1
for i=2:size(x,1)-1
  
    if x(i,j)>x(i+1,j)
        if x(i,j)>x(i-1,j)
            if x(i,j)>x(i,j-1)
                if x(i,j)>x(i,j+1)
                   y(i,j)=x(i,j);
               end
           end
       end
   end


  if x(i,j)<x(i+1,j)
        if x(i,j)<x(i-1,j)
            if x(i,j)<x(i,j-1)
               if x(i,j)<x(i,j+1)
                  v(i,j)=x(i,j);
               end
            end
        end
  end
end  
end  


 case '8'
  
for j=2:size(x,2)-1
for i=2:size(x,1)-1
  
    if x(i,j)>x(i+1,j)
      if x(i,j)>x(i+1,j+1)
        if x(i,j)>x(i+1,j-1)
          if x(i,j)>x(i,j+1)
            if x(i,j)>x(i,j-1)
              if x(i,j)>x(i-1,j)
                if x(i,j)>x(i-1,j+1)
                  if x(i,j)>x(i-1,j-1)
                   y(i,j)=x(i,j);
                  end
                end
              end
            end
          end
        end
      end
    end


    if x(i,j)<x(i+1,j)
      if x(i,j)<x(i+1,j+1)
        if x(i,j)<x(i+1,j-1)
          if x(i,j)<x(i,j+1)
            if x(i,j)<x(i,j-1)
              if x(i,j)<x(i-1,j)
                if x(i,j)<x(i-1,j+1)
                  if x(i,j)<x(i-1,j-1)
                   v(i,j)=x(i,j);
                  end
                end
              end
            end
          end
        end
      end
    end
end  
end  
case '4m'
    
    P = imregionalmax(x,4);
    y = double(P) .* double(x);
    P = imregionalmin(x,4);
    v = double(P) .* double(x);
    
    case '8m'
    
    P = imregionalmax(x,8);
    y = double(P) .* double(x);
    P = imregionalmin(x,8);
    v = double(P) .* double(x);
end

max=y;
min=v;

%*********************************************************************************************************