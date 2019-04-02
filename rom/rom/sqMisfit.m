function [p_cf_exp, Tc, TcTcT] =...
    sqMisfit(X, transType, transLimits, mesh, Tf_n_minus_mu, W_cf_n, rf2fem)
%function computing the squared difference between prediction and truth
%   X:  transformed conductivity (row vector)

%transformed conductivity to conductivity
conductivity = conductivityBackTransform(rf2fem*X', transType, transLimits);

isotropicDiffusivity = true;
if isotropicDiffusivity
    FEMout = heat2d(mesh, conductivity);
else
    %Set up conductivity tensors for each element
    D = zeros(2, 2, mesh.nEl);
    for j = 1:mesh.nEl
        D(:, :, j) =  conductivity(j)*eye(2);
    end
    %Solve coarse FEM model
    FEMout = heat2d(mesh, D);
end

Tc = FEMout.u;

p_cf_exp = (Tf_n_minus_mu - W_cf_n*Tc).^2;

if nargout > 2
    TcTcT = Tc*Tc';
end
end

