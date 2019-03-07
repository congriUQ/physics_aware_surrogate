function [sampleFun] = genBochnerSamples(lengthScale, sigma_f2, nBasisFunctions, type)
%Generate approximate Gaussian process sample functions in analytical form using Bochner's theorem

if numel(lengthScale) == 1
    lengthScale = [lengthScale lengthScale];
end

if strcmp(type, 'squaredExponential')
    %Stacked samples from W, see reference_notes
    W = mvnrnd(zeros(1, 2), diag(lengthScale.^(-2)), nBasisFunctions);
elseif strcmp(type, 'ornsteinUhlenbeck')
    %modulus
    W = trnd(2, nBasisFunctions, 1)/lengthScale(1);
    %angle
    phi = 2*pi*rand(nBasisFunctions, 1);
    W = [W.*cos(phi), W.*sin(phi)];
elseif strcmp(type, 'matern')
    %modulus - 'params' is smoothness parameter nu of Matern kernel
    %Abuse second length scale param as smoothness param
    W = trnd(lengthScale(2) + .5, nBasisFunctions, 1)/lengthScale(1);
    %angle
    phi = 2*pi*rand(nBasisFunctions, 1);
    W = [W.*cos(phi), W.*sin(phi)];
elseif strcmp(type, 'sincCov')
    W = (rand(nBasisFunctions, 1) - .5)/lengthScale(1);
    W = [W, (rand(nBasisFunctions, 1) - .5)/lengthScale(2)];
elseif strcmp(type, 'sincSqCov')
    W = triangRand(nBasisFunctions)/lengthScale(1);
    W = [W, triangRand(nBasisFunctions)/lengthScale(2)];
elseif strcmp(type, 'cos')
    Wtemp1 = randi(2, nBasisFunctions, 1);
    Wtemp2 = randi(2, nBasisFunctions, 1);
    W1(Wtemp1 == 1) = -(1/lengthScale(1));
    W1(Wtemp1 == 2) = (1/lengthScale(1));
    W2(Wtemp2 == 1) = -(1/lengthScale(2));
    W2(Wtemp2 == 2) = (1/lengthScale(2));
    W = [W1; W2]';
else
    error('Unknown covariance type')
end

%Stacked samples from b, see notes
b = 2*pi*rand(nBasisFunctions, 1);

%Draw coefficients gamma
gamma = normrnd(0, 1, 1, nBasisFunctions);

%Handle to sample function
sampleFun = @(x) sqrt((2*sigma_f2)/nBasisFunctions)*(gamma*cos(W*x + b));


end

