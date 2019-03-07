function [x] = triangRand(nSamples)
%Draws standard triangular distributed random numbers

r = rand(nSamples, 1);
x = zeros(nSamples, 1);
x(r <= .5) = -1 + sqrt(2*r(r <= .5));
x(r > .5) = 1 - sqrt(2 - 2*r(r > .5));


end

