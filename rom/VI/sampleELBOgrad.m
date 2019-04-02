function [ELBOgrad, ELBOgradErr] =...
    sampleELBOgrad(log_emp_dist, variationalDist, nSamples, varDistParams)
%Estimation of gradient of evidence lower bound (ELBO)
%log_emp_dist is a function handle to the log empirical distribution and grad

dim = numel(varDistParams.mu);
if strcmp(variationalDist, 'diagonalGauss')
    d_mu_mean = 0;
    d_sigma_mean = 0;
    d_muSq_mean = 0;
    d_sigmaSq_mean = 0;
    
    for i = 1:nSamples
        sample = normrnd(0, 1, 1, dim);
        
        %transform standard normal sample to sample of VI distribution
        variationalSample = varDistParams.mu + varDistParams.sigma.*sample;
        
        %Gradient w.r.t. dependent variable x
        [~, d_log_empirical] = log_emp_dist(variationalSample);
        d_log_empirical = d_log_empirical';
        
        %Mean gradient w.r.t. mu; d/dmu = (dX/dmu)*d/dX, X = mu + sigma*sample
        % --> d/dmu = d/dX
        d_mu_mean = (1/i)*((i - 1)*d_mu_mean + d_log_empirical);
        
        %Mean gradient w.r.t. d/d_sigma_k; d/dsigma = (dX/dsigma)*d/dX,
        %X = mu + sigma*sample --> d/dsigma = sample*(d/dX)
        %second term is due to gradient of variational dist (given analytically)
        d_sigma_mean = (1/i)*((i - 1)*d_sigma_mean +...
            (d_log_empirical.*sample + 1./varDistParams.sigma));
        
        %Might be untrue as gradient of variational distribution is missing
        %w.r.t. mu
        d_muSq_mean = (1/i)*((i - 1)*d_muSq_mean + (d_log_empirical).^2);
        %w.r.t. d/d_sigma_k
        d_sigmaSq_mean = (1/i)*((i - 1)*d_sigmaSq_mean +...
            (d_log_empirical.*sample + 1./varDistParams.sigma).^2);
    end
    
    %Transformation d/dsigma --> d/dlog(sigma^-2)
    d_logSigma_Minus2mean = -.5*(d_sigma_mean.*varDistParams.sigma);
    ELBOgrad = [d_mu_mean d_logSigma_Minus2mean];
    
    d_muErr = sqrt(abs(d_muSq_mean - d_mu_mean.^2))/sqrt(nSamples);
    %error w.r.t. d/d_sigma_k
    d_sigmaErr = sqrt(.25*(varDistParams.sigma.^2).*d_sigmaSq_mean...
        - d_logSigma_Minus2mean.^2)/sqrt(nSamples);
    ELBOgradErr = [d_muErr d_sigmaErr];
    
elseif strcmp(variationalDist, 'fullRankGauss')
    %Sigma = L*L^T
    d_mu_mean = 0;
    d_L_mean = 0;
    d_muSq_mean = 0;
    d_LSq_mean = 0;
    
    for i = 1:ELBOgradParams.nSamples
        sample = normrnd(0, 1, 1, dim);
        
        %transform standard normal sample to sample of VI distribution
        variationalSample = varDistParams.mu + sample*varDistParams.LT;
        
        %Gradient w.r.t. dependent variable x
        [~, d_log_empirical] = log_emp_dist(variationalSample);
        d_log_empirical = d_log_empirical';
        
        %Mean gradient w.r.t. mu; d/dmu = (dX/dmu)*d/dX, X = mu + sigma*sample
        % --> d/dmu = d/dX
        d_mu_mean = (1/i)*((i - 1)*d_mu_mean + d_log_empirical);
        
        %Mean gradient w.r.t. d/d_sigma_k; d/dsigma = (dX/dsigma)*d/dX,
        %X = mu + sigma*sample --> d/dsigma = sample*(d/dX)
        %second term is due to gradient of variational dist (given analytically)
        d_L_mean = (1/i)*((i - 1)*d_L_mean +...
            .5*(d_log_empirical'*sample + varDistParams.LInv)...
            + .5*(d_log_empirical'*sample + varDistParams.LInv)');
        
        %Might be untrue as gradient of variational distribution is missing
        %w.r.t. mu
        d_muSq_mean = (1/i)*((i - 1)*d_muSq_mean + (d_log_empirical).^2);
        %w.r.t. d/d_sigma_k
        d_LSq_mean = (1/i)*((i - 1)*d_LSq_mean +...
            (d_log_empirical'*sample + varDistParams.LInv).^2);
    end
    
    %We want to fix upper triangular elements to 0
    d_L_mean = tril(d_L_mean);
    d_LSq_mean = tril(d_L_mean);
    
    %Transformation d/dsigma --> d/dlog(sigma^-2)
    ELBOgrad = [d_mu_mean d_L_mean(:)'];
    
    d_muErr =...
        sqrt(abs(d_muSq_mean - d_mu_mean.^2))/sqrt(ELBOgradParams.nSamples);
    %error w.r.t. d/d_L
    d_LErr = sqrt(d_LSq_mean - d_L_mean.^2)/sqrt(ELBOgradParams.nSamples);
    ELBOgradErr = [d_muErr d_LErr(:)'];
else
    error('Unknown variational distribution')
end
end

