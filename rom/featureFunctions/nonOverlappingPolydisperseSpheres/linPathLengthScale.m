function [out] = linPathLengthScale(diskCenters, diskRadii, grd, distances)
%Gives back the parameters a and b from the theoretical lineal path model
%L(z) = a*exp(b*z)

L = zeros(grd.nCells, numel(distances));
for i = 1:numel(distances)
    L(:, i) = matrixLinealPath(diskCenters, diskRadii, grd, distances(i));
end
L = L + eps;
out = zeros(grd.nCells, 1);
for k = 1:grd.nCells
    f = fit(distances', L(k, :)', 'exp1');
    out(k) = f.b;
end