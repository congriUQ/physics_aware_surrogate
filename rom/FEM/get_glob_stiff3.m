function [K] = get_glob_stiff3(d_Ke, conductivity, nEq)
%Fast stiffness matrix assembly

K = reshape(d_Ke*sparse(conductivity), nEq, nEq);
end

