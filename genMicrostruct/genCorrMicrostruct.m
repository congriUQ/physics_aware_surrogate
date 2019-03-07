%% Script to generate list of circular exclusions distributed according to
% a warped Gaussian process
clear;
rng('shuffle');

mode = 'GP_GPR';    %engineered, GP or GP_GPR (Gaussian process also on radii)

lengthScale = .08;
lengthScale_r = .05;
sigmaGP_r = .01;
covarianceFunction = 'squaredExponential';
sigmoid_scale = 1.2;         %sigmoid warping length scale param
nBochnerSamples = 5000;
nExclusionParams = [7.8, .2];
margins = [.003, .003, .003, .003];
rParams = [-5.53, .2];
max_trials = 100;
nMeshes = 0:2;
t_max = 3600;                 %in seconds
plt = false;
resolution = 1024;
sav = true;
origin_rejection = .03; %false if no origin rejection


%% Set up save path
if strcmp(mode, 'GP')
    savepath = strcat('./data/meshSize=256/nonOverlappingDisks/margins=');
    savepath = strcat(savepath, num2str(margins(1)), '_', num2str(margins(2)),...
        '_', num2str(margins(3)), '_', num2str(margins(4)), '/N~logn/mu=', ...
        num2str(nExclusionParams(1)), '/sigma=', num2str(nExclusionParams(2)), ...
        '/x~GP/cov=', covarianceFunction, '/l=', num2str(lengthScale),...
        '/sig_scale=', num2str(sigmoid_scale), '/r~logn/mu=', num2str(rParams(1)), ...
        '/sigma=', num2str(rParams(2)));
elseif strcmp(mode, 'GP_GPR')
    savepath = strcat('./data/meshSize=256/nonOverlappingDisks/margins=');
    savepath = strcat(savepath, num2str(margins(1)), '_', num2str(margins(2)),...
        '_', num2str(margins(3)), '_', num2str(margins(4)), '/N~logn/mu=', ...
        num2str(nExclusionParams(1)), '/sigma=', num2str(nExclusionParams(2)), ...
        '/x~GP/cov=', covarianceFunction, '/l=', num2str(lengthScale),...
        '/sig_scale=', num2str(sigmoid_scale), '/r~lognGP/mu=', num2str(rParams(1)), ...
        '/sigma=', num2str(rParams(2)), '/sigmaGP_r=', num2str(sigmaGP_r),...
        '/l=', num2str(lengthScale_r));
elseif strcmp(mode, 'engineered')
    savepath = strcat('./data/meshSize=256/nonOverlappingDisks/margins=');
    savepath = strcat(savepath, num2str(margins(1)), '_', num2str(margins(2)),...
        '_', num2str(margins(3)), '_', num2str(margins(4)), '/N~logn/mu=', ...
        num2str(nExclusionParams(1)), '/sigma=', num2str(nExclusionParams(2)), ...
        '/x~engineered','/r~logn/mu=', num2str(rParams(1)), ...
        '/sigma=', num2str(rParams(2)));
elseif strcmp(mode, 'tiles')
    savepath = strcat('./data/meshSize=256/nonOverlappingDisks/margins=');
    savepath = strcat(savepath, num2str(margins(1)), '_', num2str(margins(2)),...
        '_', num2str(margins(3)), '_', num2str(margins(4)), '/N~logn/mu=', ...
        num2str(nExclusionParams(1)), '/sigma=', num2str(nExclusionParams(2)), ...
        '/x~tiles','/r~logn/mu=', num2str(rParams(1)), ...
        '/sigma=', num2str(rParams(2)));
else
    error('unknown mode');
end
if origin_rejection
    savepath = strcat(savepath, '/origin_rejection=', num2str(origin_rejection));
end
savepath

if ~exist(savepath, 'dir')
    mkdir(savepath);
end

mesh_iter = 1;
for n = nMeshes
    if(strcmp(mode, 'GP') || strcmp(mode, 'GP_GPR'))
        rejectionFun = genBochnerSamples(...
            lengthScale, 1, nBochnerSamples, covarianceFunction);
        rejectionFun = @(x) sigmf_own(rejectionFun(x), [sigmoid_scale, 0]);
    elseif strcmp(mode, 'engineered')
        rejectionFun = @(x) engineeredRejectionFun(x);
    elseif strcmp(mode, 'tiles')
        rejections = get_rejections(8);
        rejectionFun = @(x) tilesRejectionFun(x, rejections, origin_rejection);
    else
        error('unknown mode')
    end
    
    nExclusions = round(lognrnd(nExclusionParams(1), nExclusionParams(2)));
    diskCenters = zeros(nExclusions, 2);
    diskRadii = zeros(1, nExclusions);
    currentDisks = 0;
    
    t0 = tic;
    t_elapsed = 0;
    while(currentDisks < nExclusions && t_elapsed < t_max)
        diskCenter = rand(1, 2);
        if strcmp(mode, 'GP_GPR')
            radiiFun = genBochnerSamples(...
            lengthScale_r, sigmaGP_r^2, nBochnerSamples, 'squaredExponential');
            mu_r = radiiFun(diskCenter') + rParams(1);
            diskRadius = lognrnd(mu_r, rParams(2));
        else
            diskRadius = lognrnd(rParams(1), rParams(2));
        end
        
        %accept/reject
        if(rand <= rejectionFun(diskCenter'))
            %check if new circle overlaps with another circle
            
            overlap = false;
            if currentDisks     %first disk cannot overlap            
                overlap = any(((diskRadius + diskRadii(1:currentDisks))').^2 >=...
                    sum((diskCenter - diskCenters(1:currentDisks, :)).^2, 2));
                trials = 0;
                while(overlap && trials < max_trials)
                    if strcmp(mode, 'GP_GPR')
                        diskRadius = lognrnd(mu_r, rParams(2));
                    else
                        diskRadius = lognrnd(rParams(1), rParams(2));
                    end
                    overlap = any(((diskRadius + diskRadii(1:currentDisks))').^2 >=...
                        sum((diskCenter - diskCenters(1:currentDisks, :)).^2, 2));
                    trials = trials + 1;
                end
            end
            
            if ~overlap
                onBoundary = false;
                if(((diskCenter(2) - diskRadius) < margins(1) &&...
                        margins(1) >= 0) || ...
                        (((diskCenter(1) + diskRadius) > 1 - margins(2))...
                        && margins(2) >= 0) ||...
                        (((diskCenter(2) + diskRadius) > 1 - margins(3))...
                        && margins(3) >= 0) ||...
                        ((diskCenter(1) - diskRadius) < margins(4))...
                        && (margins(4) >= 0))
                    onBoundary = true;
                end
                if ~onBoundary
                    diskCenters(currentDisks + 1, :) = diskCenter;
                    diskRadii(currentDisks + 1) = diskRadius;
                    currentDisks = currentDisks + 1;
                end
            end
        end
        t_elapsed = toc(t0);
    end
    
    if currentDisks < nExclusions
        diskCenters((currentDisks + 1):end, :) = [];
        diskRadii((currentDisks + 1):end) = [];
    end
    
    %% save microstructure data
    if sav
        save(strcat(savepath, '/microstructureInformation_nomesh', num2str(n)),...
            'diskCenters', 'diskRadii');
        disp('microstructure saved.')
    end
    
    %% Plotting
    if(plt && feature('ShowFigureWindows'))      
        [img_handle, fig_handle] =...
            plotMicrostruct(diskCenters, diskRadii, resolution);
    end
    mesh_iter = mesh_iter + 1;
end


function [r] = engineeredRejectionFun(x)
    %[.0 .0]-[.25 .25] cell is empty
    if(x(1) <= .25 && x(2) <= .25)
        r = 0;
    elseif(x(1) > .25 && x(1) <= .5 && x(2) <= .25)
        r = 4*x(1) - 1;
    elseif(x(2) > .25 && x(2) <= .5 && x(1) <= .25)
        r = 4*x(2) - 1;
    elseif(x(1) > .25 && x(1) <= .5 && x(2) > .25 && x(2) <= .5)
        r = -3 + 8*x(1) + 8*x(2) - 16*x(1)*x(2);
    elseif(x(1) > .5 || x(2) > .5)
        r = 1;
    else
        r = 0;
    end

end


function [r] = tilesRejectionFun(x, rejections, origin_rejection)
    N = sqrt(numel(rejections));
    [xx, yy] = meshgrid(linspace(0, 1, N + 1));
    xxl = xx(1:N, 1:N);
    yyl = yy(1:N, 1:N);
    xxu = xx(1:N, 2:(N + 1));
    yyu = yy(2:(N + 1), 1:N);
    l = [xxl(:), yyl(:)];
    u = [xxu(:), yyu(:)];

    cell = find((x(1) >= l(:, 1)).*(x(2) >= l(:, 2)).*(x(1) <...
        u(:, 1)).*(x(2) < u(:, 2)));
    r = rejections(cell);
    if norm(x) < origin_rejection
        r = 0;
    end
end


function [r] = get_rejections(N)
    %[.0 .0]-[.25 .25] cell is empty
    %N is the linear cell number i.e. 8 for an 8x8 mesh
    [xx, yy] = meshgrid(linspace(0, 1, N + 1));
    xxu = xx(1:N, 2:(N + 1));
    yyu = yy(2:(N + 1), 1:N);
    u = [xxu(:), yyu(:)];
    
    r = .3*ones(64, 1);
    
    for i = 1:64
        if(u(i, 1) <= .5 && u(i, 2) <= .5)
            %lower left square
            r(i) = rand;
        end
    end
end


function r = sigmf_own(x, params)
    r = 1./(1 + exp(-params(1)*(x - params(2))));
end













