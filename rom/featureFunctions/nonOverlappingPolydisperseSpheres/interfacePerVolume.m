function [relativeInterfaceArea] = interfacePerVolume(...
    diskCenters, diskRadii, gridRF)
%Computes the volume fraction of microstructures with disk exclusions
%   diskCenters:         clear
%   diskRadii:           clear
%   gridVectorX/Y:       specification of macro-cell edge lengths


A0 = zeros(gridRF.nCells, 1);
interface_area = A0;
circle_circumferences = 2*pi*diskRadii;
n = 1;
for cll = gridRF.cells
    if isvalid(cll{1})
        A0(n) = cll{1}.surface;
        circles_in_n = cll{1}.inside(diskCenters);
        interface_area(n) = sum(circle_circumferences(circles_in_n));
        n = n + 1;
    end
end

relativeInterfaceArea = interface_area./A0;




