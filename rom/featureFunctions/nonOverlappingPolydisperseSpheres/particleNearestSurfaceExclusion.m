function [e_p, h_p, opt_d, opt_h_p] = particleNearestSurfaceExclusion(...
    diskCenters, diskRadii, gridRF, distance, ref_distance)
%Mean chord length for non-overlapping polydisp. spheres
%according to Torquato 6.46, 6.47
% ref_distance should be smaller than distance

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
meanRadii(~isfinite(meanRadii)) = 0;
meanSqRadii(~isfinite(meanSqRadii)) = 0;

porefrac = A./A0;
porefrac(porefrac <= 0) = eps;  %can happen if circles lie on boundary

exclfrac = 1 - porefrac;

S = (meanRadii.^2)./(meanSqRadii + eps);
a_0 = (1 + exclfrac.*(S - 1))./(porefrac.^2);
a_1 = 1./porefrac;

x = distance./(2*meanRadii);
x(~isfinite(x)) = 0;
X = ref_distance./(2*meanRadii);
X(~isfinite(X)) = 0;
e_p = exp(-4*exclfrac.*S.*(a_0.*(x.^2 - X.^2) + a_1.*(x - X)));
if any(~isfinite(e_p))
    %warning('Non-finite value in voidNearestSurfaceExclusion')
    e_p(~isfinite(e_p)) = 1;
end
if nargout > 1
    h_p = 2*((exclfrac.*S)./(meanRadii)).*(2*a_0.*x + a_1).*e_p;
end
if any(~isfinite(h_p))
    %warning('Non-finite value in voidNearestSurfaceExclusion')
    h_p(~isfinite(h_p)) = 0;
end

if nargout > 2
    %maximum of h_p, ignore const. prefactors!
    % + eps in denominator for numerical stability
    xx = @(d) d./(2*meanRadii + eps);
    e_p_fun =...
        @(d) exp(-4*exclfrac.*S.*(a_0.*(xx(d).^2 - X.^2) + a_1.*(xx(d) - X)));
    neg_h_p = @(d) -(2*a_0.*xx(d)+a_1).*e_p_fun(d);

    opt_d = zeros(gridRF.nCells, 1);
    opt_h_p = opt_d;
    pick = opt_d;
    for k = 1:gridRF.nCells
        pick = 0*pick;
        pick(k) = 1;
        fun = @(d) pick'*neg_h_p(d);
        [opt_d(k), opt_h_p(k)] = fminsearch(fun, .01);
    end
    exception_cells = opt_d > sqrt(A0);
    opt_d(exception_cells) = sqrt(A0(exception_cells));
    opt_h_p(exception_cells) = 0;
end


