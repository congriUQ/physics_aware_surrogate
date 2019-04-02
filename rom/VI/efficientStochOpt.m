function [varDistParams, x] = efficientStochOpt(...
    x, log_emp_dist, variationalDist, stepWidth, dim, maxCompTime)
%Memory efficient stocastic optimization for parfor loop
%Perform stochastic maximization step

debug = false;   %debug mode

updateRule = 'amsgrad';
% beta1 = .7;                     %the higher, the more important is momentum
% beta2 = .8;                    %curvature parameter
beta1 = .9;
beta2 = .9999995;
epsilon = 1e-6;                  %curvature stabilization parameter

stepOffset = 200000;                %Robbins-Monro step offset
maxIterations = Inf;
if nargin < 6
    maxCompTime = 10;
end
nSamplesStart = 1;                  %gradient samples per iteration
nSamplesEnd = 1;
nIncr = (nSamplesEnd - nSamplesStart)/maxCompTime;
nSamples = nSamplesStart;

converged = false;
steps = 0;

stepWidth_stepOffset = stepWidth*stepOffset;
if strcmp(variationalDist, 'diagonalGauss')
    varDistParams.mu = x(1:dim);
    varDistParams.sigma = exp(-.5*x((dim + 1):end));
elseif strcmp(variationalDist, 'fullRankGauss')
    varDistParams.mu = x(1:dim);
    %lower triangular cholesky factor L
    varDistParams.L = reshape(x((dim + 1):end), dim, dim);
    varDistParams.LT = varDistParams.L';
    varDistParams.LInv = inv(varDistParams.L);
    varDistParams.Sigma = varDistParams.L*varDistParams.LT;
else
    error('Unknown variational distribution')
end

if debug
    disp('starting stochastic optimization')
    x_0 = x
    f = figure;
    f.Units = 'normalized';
    f.OuterPosition = [0 0 1 1];
    clf(f)
    sb1 = subplot(2, 2, 1, 'Parent', f);
    hold on;
    title('$\mu$');
    sb2 = subplot(2, 2, 2, 'Parent', f);
    sb2.YScale = 'log';
    hold on;
    title('$\sigma$')
    sb3 = subplot(2, 2, 3, 'Parent', f);
    sb3.YScale = 'log';
    hold on;
    title('norm momentum');
    sb4 = subplot(2, 2, 4, 'Parent', f);
    sb4.YScale = 'log';
    hold on;
    title('norm gradient');
    
    for d = 1:dim
        p_mu(d) = animatedline('color', 'b', 'Parent', sb1);
        p_sigma(d) = animatedline('color', 'r', 'Parent', sb2);
    end
    p_momentum = animatedline('Parent', sb3);
    p_grad = animatedline('Parent', sb4);
    
    t_plot = tic;
end

cmpt = tic;
while ~converged
    
    [gradient] =...
        sampleELBOgrad(log_emp_dist, variationalDist, nSamples, varDistParams);
    

    if strcmp(updateRule, 'adam')
        
        if steps == 0
            %careful first iteration
            momentum = 1e-6*gradient;
            uncenteredXVariance = gradient.^2;
        else
            momentum = beta1*momentum + (1 - beta1)*gradient;
        end
        uncenteredXVariance = beta2*uncenteredXVariance...
            + (1 - beta2)*gradient.^2;
        
        %Optimization update
        x = x + (stepWidth_stepOffset/(stepOffset + steps)).*...
            (1./(sqrt(uncenteredXVariance) + epsilon)).*momentum;
        
    elseif strcmp(updateRule, 'amsgrad')
        
        if steps == 0
            %careful first iteration
            momentum = 1e-6*gradient;
            uncenteredXVariance = gradient.^2;
            uncenteredXVariance_max = uncenteredXVariance;
        else
            momentum = beta1*momentum + (1 - beta1)*gradient;
        end
        uncenteredXVariance = beta2*uncenteredXVariance...
            + (1 - beta2)*gradient.^2;
        uncenteredXVariance_max(uncenteredXVariance_max<uncenteredXVariance)...
            =uncenteredXVariance(uncenteredXVariance_max < uncenteredXVariance);
        
        %Optimization update
        x = x + (stepWidth_stepOffset/(stepOffset + steps)).*...
            (1./(sqrt(uncenteredXVariance_max) + epsilon)).*momentum;
        
    elseif strcmp(updateRule, 'robbinsMonro')
        delta = ((stepWidth_stepOffset)/(stepOffset + steps)).*gradient;
        nDelta = norm(delta);
        stabilityFactor = 2;
        if(nDelta > stabilityFactor*norm(x))
            delta = (stabilityFactor/nDelta)*delta;
        end
        x = x + delta;
    else
        error('Unknown update heuristic for stochastic optimization')
    end
    steps = steps + 1;
    
    if strcmp(variationalDist, 'diagonalGauss')
        varDistParams.mu = x(1:dim);
        varDistParams.sigma = exp(-.5*x((dim + 1):end));
    elseif strcmp(variationalDist, 'fullRankGauss')
        varDistParams.mu = x(1:dim);
        %lower triangular cholesky factor L
        varDistParams.L = reshape(x((dim + 1):end), dim, dim);
        varDistParams.LT = varDistParams.L';
        varDistParams.LInv = inv(varDistParams.L);
        varDistParams.Sigma = varDistParams.L*varDistParams.LT;
    else
        error('Unknown variational distribution')
    end
    
    if debug
        t = toc(t_plot);
        compTime = toc(cmpt);
        if(t > 2)
            addpoints(p_grad, compTime, norm(gradient));
            addpoints(p_momentum, compTime, norm(momentum));
            for d = 1:dim
                addpoints(p_mu(d), compTime, varDistParams.mu(d));
                addpoints(p_sigma(d), compTime, varDistParams.sigma(d));
            end
            axis(sb1, 'tight'), axis(sb1, 'fill');
            axis(sb2, 'tight'), axis(sb2, 'fill');
            axis(sb3, 'tight'), axis(sb3, 'fill');
            axis(sb4, 'tight'), axis(sb4, 'fill');
            drawnow
            t_plot = tic;
        end
    end
    
    compTime = toc(cmpt);
    nSamples = ceil(nIncr*compTime + nSamplesStart);
    if steps > maxIterations
        converged = true;
        disp('Converged because max number of iterations exceeded')
    elseif compTime > maxCompTime
        converged = true;
        disp('Converged because max computation time exceeded')
        if debug
            steps
        end
    end
end

if strcmp(variationalDist, 'diagonalGauss')
    varDistParams.XSqMean = varDistParams.sigma.^2 + varDistParams.mu.^2;
elseif strcmp(variationalDist, 'fullRankGauss')
    varDistParams.XSqMean =...
        diag(varDistParams.Sigma + varDistParams.mu'*varDistParams.mu);
else
    error('Unknown variational distribution')
end


