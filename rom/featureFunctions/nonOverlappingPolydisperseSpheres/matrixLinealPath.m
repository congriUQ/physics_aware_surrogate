function [L] = matrixLinealPath(diskCenters, diskRadii,...
    gridRF, distance)
%Gives the lineal path function for the matrix phase according to the
%approximation in Torquato eq. 6.37


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
porefrac(porefrac <= 0) = eps;  %can happen if circles lie on boundary


L = porefrac.*...
    exp(-(2*distance*(1 - porefrac).*meanRadii)./(pi*porefrac.*meanSqRadii));
L(isnan(L)) = 1; %Happens if there are no exclusions in macro-cell



end

