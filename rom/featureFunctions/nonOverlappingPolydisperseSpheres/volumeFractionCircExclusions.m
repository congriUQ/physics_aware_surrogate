function [porefrac] = volumeFractionCircExclusions(...
    diskCenters, diskRadii, gridRF)
%Computes the pore fraction of microstructures with disk exclusions
%pores == voids where fluid can flow
%   diskCenters:         clear
%   diskRadii:           clear
%   gridRF:              random field mesh

A = zeros(gridRF.nCells, 1);
A0 = A;
circle_surfaces = pi*diskRadii.^2;
n = 1;
for cll = gridRF.cells
    if isvalid(cll{1})  %check if cell has been deleted during splitting
        A(n) = cll{1}.surface;
        A0(n) = A(n);
        circles_in_n = cll{1}.inside(diskCenters);
        A(n) = A0(n) - sum(circle_surfaces(circles_in_n));
        n = n + 1;
    end
end

porefrac = A./A0;
porefrac(porefrac <= 0) = eps;  %can happen if circles lie on boundary

end

