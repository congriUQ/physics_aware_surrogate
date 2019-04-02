function [cld] = chordLengthDensity(diskCenters, diskRadii, gridRF, distance)
%Mean chord length for non-overlapping polydisp. spheres
%according to Torquato 6.39

meanRadii = zeros(gridRF.nCells, 1);
meanSqRadii = zeros(gridRF.nCells, 1);
A = zeros(gridRF.nCells, 1);
A0 = A;

n = 1;
for cll = gridRF.cells
    if isvalid(cll{1})  %check if cell has been deleted during splitting
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
exclfrac = 1 - porefrac;


lc = .5*pi*(meanSqRadii./meanRadii).*(porefrac./exclfrac);
if any(~isfinite(lc))
    %warning('Setting mean chord length of cells with no inclusion to 0.')
    lc(isnan(lc)) = sqrt(A(isnan(lc))); %set to ~ cell edge length
end

cld = (1./lc).*exp(-(distance./lc));