function [F] = get_glob_force(mesh, conductivity)
%Assemble global force vector

F = mesh.F_natural;
% for performance
anyEssentialNodeInElement = any(mesh.essentialNodeInElement');
for e = 1:mesh.nEl
%     f = get_loc_force(e, domain, k);
    %Contribution due to essential boundaries
    %local stiffness matrix k

    if anyEssentialNodeInElement(e)
        Tb = zeros(4, 1);
        for i = 1:4
            if mesh.essentialNodeInElement(e, i)
                Tb(i) = mesh.essentialTemperatures(mesh.globalNodeNumber(e, i));
            end
        end

        for ln = 1:4
            eqn = mesh.lm(e, ln);
            if(eqn ~= 0)
                fT = conductivity(e)*mesh.d_loc_stiff(:, :, e)*Tb;
                F(eqn) = F(eqn) - fT(ln);
            end
        end
    end
end

end

