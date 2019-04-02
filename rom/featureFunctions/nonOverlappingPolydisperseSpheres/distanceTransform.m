function [out] = distanceTransform(...
    lambdaMat, distMeasure, phaseInversion, meanVarMaxMin)
%Uses built-in Matlab 'bwdist' to return mean/var/max/min of pixel distances to
%next phase
%   lambdaMat:          2-dim conductivity image
%   conductivities:     loCond in first, upCond in second entry
%   hilo:               property for high or low phase bubbles?
%   distMeasure: 'euclidean', 'chessboard', 'cityblock', 'quasi-euclidean'
%See matlab reference for bwdist

if phaseInversion
    lambdaMat = ~lambdaMat;
end

dist = bwdist(lambdaMat, distMeasure);

%Catch infinities: Take maximum possible distance
if(any(any(isinf(dist))))
%     warning(...
%'Infinity in distance transformation. Setting to maximum possible distance')
    if strcmp(distMeasure, 'cityblock')
        dist(isinf(dist)) = size(dist, 1) + size(dist, 2);
    elseif strcmp(distMeasure, 'chessboard')
        dist(isinf(dist)) = max([size(dist, 1), size(dist, 2)]);
    else
        %Treat euclidean and quasi-euclidean equally. This is actually wrong
        %for quasi-euclidean
        dist(isinf(dist)) = norm(size(dist));
    end
end

if strcmp(meanVarMaxMin, 'mean')
    out = mean(mean(dist));
elseif strcmp(meanVarMaxMin, 'var')
    out = var(dist(:));
elseif strcmp(meanVarMaxMin, 'max')
    out = max(max(dist));
elseif strcmp(meanVarMaxMin, 'min')
    warning('Min dist is usually 0 for every cell. Feature should not be used.')
    out = min(min(dist));
else
    error('Mean or variance of lambda bubble property?')
end


end

