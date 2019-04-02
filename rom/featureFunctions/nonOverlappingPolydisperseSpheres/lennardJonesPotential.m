function [V] = lennardJonesPotential(distances, d)
%Lennard-Jones potential energy for spherical exclusions

V = sum((d./distances).^12 - (d./distances).^6);

end