function [V] = squareWellPotential(distances, d)
%Computes square well potential for microstructure exclusions
%   d:  potential width

V = -sum(distances(distances < d));
end