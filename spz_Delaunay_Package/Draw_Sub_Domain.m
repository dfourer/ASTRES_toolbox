function [mask, borders] = Draw_Sub_Domain(triZ, xZ, yZ, Seg_Index_Vector)



x_min = min(min(xZ(triZ(Seg_Index_Vector, :))));
x_max = max(max(xZ(triZ(Seg_Index_Vector, :))));
y_min = min(min(yZ(triZ(Seg_Index_Vector, :))));
y_max = max(max(yZ(triZ(Seg_Index_Vector, :))));

% fprintf('\n\tx_min: %d, \tx_max: %d,\ty_min: %d,\ty_max %d,', x_min, x_max, y_min, y_max);
borders = [x_min, x_max, y_min, y_max];
NAx = x_max - x_min +1;
NAy = y_max - y_min +1;

mask = zeros(NAx, NAy);

for this_Triangle_Index = Seg_Index_Vector
    this_Triangle = triZ(this_Triangle_Index, :);
    x = [xZ(this_Triangle); xZ(this_Triangle(1))] - x_min + 1;
    y = [yZ(this_Triangle); yZ(this_Triangle(1))] - y_min + 1;
    this_mask = poly2mask(y, x, NAx, NAy);
    mask = mask + this_mask;
end
