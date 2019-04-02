function [s] = specificSurface(diskCenters, diskRadii, gridRF)
%Specific surface for non-overlapping polydisperse spheres
%See Torquato 6.33

meanRadii = zeros(gridRF.nCells, 1);
meanSqRadii = zeros(gridRF.nCells, 1);
A = zeros(gridRF.nCells, 1);
A0 = A;

n = 1;
for cll = gridRF.cells
    if isvalid(cll{1})
        radii_in_n = diskRadii(cll{1}.inside(diskCenters));
        meanRadii(n) = mean(radii_in_n);
        meanSqRadii(n) = mean(radii_in_n.^2);
        A0(n) = cll{1}.surface;
        A(n) = cll{1}.surface - pi*sum(radii_in_n.^2);
        n = n + 1;
    end
end

porefrac = A./A0;
s = 2*(1 - porefrac).*(meanRadii./meanSqRadii);
s(isnan(s)) = 0; %this occurs if macro-cell has no inclusions

end

