function [lc] = meanChordLength(diskCenters, diskRadii, gridRF)
%Mean chord length for non-overlapping polydisp. spheres
%according to Torquato 6.40

meanRadii = zeros(gridRF.nCells, 1);
meanSqRadii = zeros(gridRF.nCells, 1);
edg_max = meanRadii;
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
        for edg = 1:numel(cll{1}.edges)
            if isvalid(cll{1}.edges{edg})
                if cll{1}.edges{edg}.length > edg_max(n)
                    edg_max(n) = cll{1}.edges{edg}.length;
                end
            else
                warning('Cell with deleted edge');
            end
        end
        n = n + 1;
    end
end

porefrac = A./A0;
porefrac(porefrac <= 0) = eps;  %can happen if circles lie on boundary
exclfrac = 1 - porefrac;

lc = .5*pi*(meanSqRadii./meanRadii).*(porefrac./exclfrac);
%Set maximum chord length to maximum edge length
lc(lc > edg_max) = edg_max(lc > edg_max);
%can happen in macro-cells with no exclusions
lc(~isfinite(lc)) = edg_max(~isfinite(lc));

