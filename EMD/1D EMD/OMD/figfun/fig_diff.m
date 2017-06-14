function [a f er] = fig_diff(ares,fres)
% Computes the function ac = h(f) which is the curve separating the domains
%   - D1 : s and s'' have the same number of extrema.
%   - D2 : s and s'''' have the same number of extrema.
%   - D3 : not (D1 or D2)
% ares, fres : horizontal and vertical resolution (default : 0.01 and
% 0.005).

if nargin<2
    ares = 0.005;
    fres = 0.01;
end


T = 8*1024;
seuil= 4;

a = 10.^(-2:ares:2);
f  = 0.01:fres:1;
t = linspace(0,T,128*T);
intval = 2*T:128*T;

er = zeros(length(a),length(f));
x1 = cos(2*pi*t);

for i =1:length(a)
 for j = 1:length(f),
     
  s = x1+a(i)*cos(2*pi*f(j)*t);
  
  sder2 = madiff(s,t,2);
  sder4 = madiff(s,t,4);
  
  
  
  s = s(intval);
  sder2 = sder2(intval);
  sder4 = sder4(intval);
  
  
  % Par rapport aux extrema de s, s'' ou s4
  [indmin, indmax, indzer] = extr(s,t);
  ls0 = length(indmin) + length(indmax);
        
  [indmin, indmax, indzer] = extr(sder2);
  ls2 = length(indmin) + length(indmax);
  
  [indmin, indmax, indzer] = extr(sder4);
  ls4 = length(indmin) + length(indmax);
  
  if ls2<=ls0+seuil
      % On choisit s
      er(i,j)=0;
  elseif ls4 <=ls2+seuil
      % On choisit s2
      er(i,j)=1;
  else
      % On choisit s4
      er(i,j)=2;
  end
  
  if max(ls0,max(ls2,ls4))>=length(s)/5
      error('Pb de discretisation dans le calculm de la derivee');
  end
        
  
 end
end

er=er';

figure();
imagesc(log(a)/log(10),(f),er);
set(gca,'YDir','normal');
colormap('gray');
xlabel('a');
ylabel('f');
hold on;
plot((-2:0.1:2),(10.^(-2:0.1:2)).^(-1/1),'color',0.3*ones(1,3),'LineWidth',2);
plot((-2:0.1:2),(10.^(-2:0.1:2)).^(-1/3),'color',0.5*ones(1,3),'LineWidth',2);
plot((-2:0.1:2),(10.^(-2:0.1:2)).^(-1/5),'color',0.7*ones(1,3),'LineWidth',2);
legend('af=1','af^3=1','af^5=1');
set(gca,'xtick',-2:2,'xticklabel',{'10^{-2}','10^{-1}','1','10^{1}','10^{2}',});


end
