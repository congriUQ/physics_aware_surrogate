function [k] = kozenyCarman(diskCenters, diskRadii, gridRF)
%See e.g. eq. 3 in Microscale permeability predictions of porous fibrous media

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


%the +1 in the denominator is for regularity for high pore fractions
k = (porefrac.^3)./(5*(s.^2).*(1 - porefrac).^2 + 1);
k(~isfinite(k)) = 1;

end