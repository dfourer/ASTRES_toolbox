function [ as ] = a_axis(M,  as_range, method, i)
% [ mm ] = a_axis(M, i)
% Compute the scale values in range [as_range(1), as_range(2)]
% 
% INPUT:
% M       : number of bins to compute
% is_freq : (optional) select how scale axis is sampled 
%          0: (default) logarithmic scale (better for signal reconstruction)
%          1: 1/f scale, with f sampled linearly (better to display using frequencies)
%             as_range is assumed to contain [w_begin/w0  w_end/w0]
%             where w_begin (resp. w_end) is equal to 2*pi*f_begin (resp. 2*pi*f_end)
%          2: linear scale (better to display using scales)
%
% i       : (optional) selected bins in
%
%
% OUTPUT:
% mm     : the computed frequency axis 
%
% Author: D.Fourer
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. Signal Processing. 2015]

if ~exist('method', 'var')
 method = 0; 
end


if method == 0        %% log scale (better for signal reconstruction)
 as = logspace(log10(as_range(1)), log10(as_range(2)), M);
   
elseif method == 1    %% 1/f scale (better for display using frequencies)
 as = 1 ./ linspace(as_range(1), as_range(2), M);
 
elseif method == 2    %% linear scale (better for display using scales)
 as = linspace(as_range(1), as_range(2), M); 
 
end
 
if exist('i', 'var')
  as = as(i);
end