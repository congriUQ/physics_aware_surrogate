function [conductivity] = conductivityBackTransform(x, transType, transLimits)
%Backtransformation from X to conductivity lambda

if strcmp(transType, 'logit')
    %Logistic sigmoid transformation
    conductivity =...
        (transLimits(2) - transLimits(1))./(1 + exp(-x)) + transLimits(1);
elseif strcmp(transType, 'log_cholesky')
    %conductivity is a 2x2 tensor now, store in 3-dim array
    N = size(x, 2);
    if(~(size(x, 1) == 3))
        error('X must contain 3 parameters for each conductivity tensor')
    end
    conductivity = zeros(2, 2, N);
    for i = 1:N
        %cholesky matrix
        L = [exp(x(1, i)), x(2, i); 0, exp(x(3, i))];
        conductivity(:, :, i) = L'*L;
    end
elseif strcmp(transType, 'log')
    conductivity = exp(x);
    if any(any(conductivity < transLimits(1)))
        %lower bound on conductivity for stability
        conductivity(conductivity < transLimits(1)) = transLimits(1);
%         warning('Effective conductivity is below lower bound')
    end
    if any(any(conductivity > transLimits(2)))
        %upper bound on conductivity for stability
        conductivity(conductivity > transLimits(2)) = transLimits(2);
%         warning('Effective conductivity is above upper bound')
    end
elseif strcmp(transType, 'log_lower_bound')
    %log transformation, where the eff. cond. lower bound is non-zero 
    conductivity = exp(x) + transLimits(1);
    if any(any(conductivity > transLimits(2)))
        %upper bound on conductivity for stability
        conductivity(conductivity > transLimits(2)) = transLimits(2);
        warning('Effective conductivity is above upper bound')
    end
elseif strcmp(transType, 'square')
    conductivity = x.^2;
    if any(any(conductivity < transLimits(1)))
        %lower bound on conductivity for stability
        conductivity(conductivity < transLimits(1)) = transLimits(1);
%         warning('Effective conductivity is below lower bound')
    end
    if any(any(conductivity > transLimits(2)))
        %upper bound on conductivity for stability
        conductivity(conductivity > transLimits(2)) = transLimits(2);
%         warning('Effective conductivity is above upper bound')
    end
else
    error('unknown conductivity transformation')
end

if(any(any(~isfinite(conductivity))))
    warning('Non-finite conductivity, setting it to 1e-3.')
    conductivity(~isfinite(conductivity)) = 1e-3;
end
end

