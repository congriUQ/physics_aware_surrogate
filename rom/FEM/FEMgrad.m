function [d_r] = FEMgrad(u, mesh)
%Compute derivatives of FEM equation system r = K*Y - F w.r.t. Lambda_e
%ONLY VALID FOR ISOTROPIC HEAT CONDUCTIVITY MATRIX D!!!
% u:    PDE response vector
% mesh: mesh object

% (d/d Lambda_e) k^(e) = (1/Lambda_e) k^(e)     as k^(e) linear in Lambda_e
d_r = reshape(mesh.d_glob_stiff*sparse(u), mesh.nEq, mesh.nEl)'...
    - mesh.d_glob_force;

end

