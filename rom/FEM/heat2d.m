function [Out] = heat2d(mesh, conductivity)
%2D heat conduction main function
%Gives back temperature on point x

%Global stiffness matrix
% Out.globalStiffness = get_glob_stiff2(mesh, Out.diffusionStiffness);
Out.globalStiffness = ...
    get_glob_stiff3(mesh.d_glob_stiff_assemble, conductivity, mesh.nEq);

%Global force vector
% Out.globalForce = get_glob_force(mesh, conductivity);

%Finally solving the equation system
Out.naturalTemperatures = Out.globalStiffness\get_glob_force(mesh,conductivity);

%Temperature field
%Out.u = zeros(mesh.nNodes, 1);
Out.u(mesh.id, 1) = Out.naturalTemperatures;
Out.u(mesh.essentialNodes, 1) = mesh.essentialTemperatures(mesh.essentialNodes);

% for i = 1:mesh.nNodes
%     if(any(i == mesh.essentialNodes))
%         %node i is essential
%         Out.Tf(i) = mesh.essentialTemperatures(i);
%     end
% end

end


