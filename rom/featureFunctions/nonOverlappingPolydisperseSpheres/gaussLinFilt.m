function [X] = gaussLinFilt(lambda, muGaussFilt, sigmaGaussFiltFactor)
%Gaussian shaped linear filter centered at the element center

if(nargin < 2 || any(isnan(muGaussFilt)))
    muGaussFilt = [(size(lambda, 1) + 1)/2 (size(lambda, 2) + 1)/2];
end
if nargin < 3
    sigmaGaussFiltFactor = 10;
end
sigmaGaussFilt = sigmaGaussFiltFactor*[size(lambda, 1) size(lambda, 2)];

[x, y] = meshgrid(1:size(lambda, 1), 1:size(lambda, 2));
xy = [x(:) y(:)];
w = mvnpdf(xy, muGaussFilt, sigmaGaussFilt);
w = w/sum(w);
w = reshape(w, size(lambda, 1), size(lambda, 2));

debug = false;
if debug
    figure
    subplot(1,2,1)
    imagesc(w)
    drawnow
    pause
end


X = sum(sum(w.*lambda));
end

