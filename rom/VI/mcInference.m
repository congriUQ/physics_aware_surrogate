function [E, E2, E3] = mcInference(functionHandle, variationalDist, varDistParams)
%Samples expected value of function in function handle under variational dist.
inferenceSamples = 2000;
if strcmp(variationalDist, 'diagonalGauss')
    %individual samples are rows of 'samples'
    %sigma is std
    samples = mvnrnd(varDistParams.mu, varDistParams.sigma, inferenceSamples);
    E = 0;
    E2 = 0;
    E3 = 0;
    for i = 1:inferenceSamples
        if nargout > 1
            [p_cf_exp, Tc, TcTcT] = functionHandle(samples(i, :));
            E2 = (1/i)*((i - 1)*E2 + Tc);
            E3 = (1/i)*((i - 1)*E3 + TcTcT);
        else
            p_cf_exp = functionHandle(samples(i, :));
        end
        E = (1/i)*((i - 1)*E + p_cf_exp);
    end
elseif strcmp(variationalDist, 'fullRankGauss')
    %individual samples are rows of 'samples'
    %Sigma is covariance
    samples = mvnrnd(varDistParams.mu, varDistParams.Sigma, inferenceSamples);
    E = 0;
    for i = 1:inferenceSamples
        E = (1/i)*((i - 1)*E + functionHandle(samples(i, :)));
    end
else
    error('Unknown variational distribution')
end
end

