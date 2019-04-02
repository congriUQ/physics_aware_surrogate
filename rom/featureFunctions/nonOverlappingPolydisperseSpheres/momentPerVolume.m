function [relativeInterfaceArea] = momentPerVolume(...
    diskCenters, diskRadii, gridRF, moment)
%Computes the sum of moments of disk radii divided by cell area
%   diskCenters:         clear
%   diskRadii:           clear
%   gridVectorX/Y:       specification of macro-cell edge lengths


A0 = zeros(gridRF.nCells, 1);
sum_radii_moments = A0;
radii_samples = diskRadii.^moment;
n = 1;
for cll = gridRF.cells
    if isvalid(cll{1})
        A0(n) = cll{1}.surface;
        circles_in_n = cll{1}.inside(diskCenters);
        sum_radii_moments(n) = sum(radii_samples(circles_in_n));
        n = n + 1;
    end
end

relativeInterfaceArea = sum_radii_moments./A0;

