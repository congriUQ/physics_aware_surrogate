classdef ModelParams < matlab.mixin.Copyable
    %Initialize, update, ... the ROM model params
    
    properties
        %% Model parameters (or expected values thereof)
        %p_c
        theta_c
        theta_cSq
        Sigma_c
        %posterior variance of theta_c, given a prior model
        Sigma_theta_c
        %FEM grid of coarse Darcy emulator
        coarseGridX = (1/16)*ones(1, 16)
        coarseGridY = (1/16)*ones(1, 16)
        %grid of random field
        gridRF
        splitted_cells
        
        %Matrix summing up all values of a fine-scale vector belonging to
        %a certain macro-cell
        sum_in_macrocell
        
        %Recorded elbo
        elbo
        cell_score
        cell_score_full   %includes also terms of p_cf to elbo cell score
        active_cells
        active_cells_S
        sigma_cf_score    %refine where sqrt(S) of p_cf is maximal
        inv_sigma_cf_score%refine where inv(sqrt(S)) is maximal
        
        %p_cf
        W_cf
        sigma_cf
        fineGridX = (1/128)*ones(1, 128)
        fineGridY = (1/128)*ones(1, 128)
        interpolationMode = 'linear'   %'nearest', 'linear' or 'natural'
        smoothingParameter = []
        boundarySmoothingPixels = -1   %only smooths boundary if positive
        
        %Surrogate FEM mesh
        coarseMesh
        %Mesh object identical to FEM mesh for computation of map2fine
        %(is ugly and should be merged with coarseMesh object)
        coarseMesh_geometry
        
        %Mapping from random field to FEM discretization
        rf2fem
        %cell index dictionary
        cell_dictionary
        
        %Transformation options of diffusivity parameter
        diffTransform = 'log'
        diffLimits = [1e-10, 1e8];
        
        %% Model hyperparameters
        mode = 'local'  %separate theta_c's per macro-cell
        %VRVM, sharedVRVM or adaptiveGaussian
        prior_theta_c = 'sharedVRVM'
        diag_theta_c = true     %diagonal Gaussian for Q(theta_c)
        gamma   %Gaussian precision of prior on theta_c
        VRVM_a = eps
        VRVM_b = eps
        VRVM_c = eps
        VRVM_d = eps
        VRVM_e = eps
        VRVM_f = eps
        %iterations with fixed q(lambda_c)
        %should be more in the end because q(lambda_c) is not changing
        %but the gamma's and other params are still changing
        VRVM_iter = [Inf]
        VRVM_time = [1*ones(1, 10), 2*ones(1, 10), ...
            5]
        
        %current parameters of variational distributions
        a
        b
        c
        d
        e
        f
        
        %% Parameters of variational distributions
        varDistParamsVec
        variational_mu
        variational_sigma
        
        %% Parameters to rescale features
        normalization = 'rescale'
        featureFunctionMean
        featureFunctionSqMean
        featureFunctionMin
        featureFunctionMax
        
        %% Training parameters
        EM_iter       = 0    %current total number of iterations
        EM_iter_split = 0    %current number of iterations after last split
        epoch = 0            %current epoch
        max_EM_epochs = 200   %maximum number of epochs
        
        %% Settings
        computeElbo = true
        
        %% plots, initialized as animated lines
        plt_params = true    %plot parameter monitoring plots?
        pParams
        pElbo
        pCellScores
    end
    
    methods
        function self = ModelParams(u_bc, p_bc)
            %Constructor; Set up all params that are unchanged during training
            %   u_bc:       boundary velocity field
            %   p_bc:       boundary pressure field
            
            %grid of random field
            self.gridRF = RectangularMesh((1/4)*ones(1, 4));
            self.cell_dictionary = 1:self.gridRF.nCells;
            
            %% Initialize coarse mesh object
            %Coarse mesh object
            self.coarseMesh = MeshFEM(self.coarseGridX, self.coarseGridY);
            self.coarseMesh.compute_grad = true;
            self.coarseMesh_geometry = RectangularMesh(self.coarseGridX,...
                self.coarseGridY);
            
            %% Set up coarse model bc's
            %Convert flow bc string to handle functions
            u_x_temp = strrep(u_bc{1}(5:end), 'x[1]', '*y');
            %This is only valid for unit square domain!!
            u_x_temp_le = strrep(u_x_temp, 'x[0]', '0');
            u_x_temp_r = strrep(u_x_temp, 'x[0]', '1');
            
            u_y_temp = strrep(u_bc{2}(5:end), 'x[0]', '*x');
            u_y_temp_lo = strrep(u_y_temp, 'x[1]', '0');
            u_y_temp_u = strrep(u_y_temp, 'x[1]', '1');
            u_bc_handle{1} = str2func(strcat('@(x)', '-(', u_y_temp_lo, ')'));
            u_bc_handle{2} = str2func(strcat('@(y)', u_x_temp_r));
            u_bc_handle{3} = str2func(strcat('@(x)', u_y_temp_u));
            u_bc_handle{4} = str2func(strcat('@(y)', '-(', u_x_temp_le, ')'));
            
            p_bc_handle = str2func(strcat('@(x)', p_bc));
            
            nX = length(self.coarseGridX);
            nY = length(self.coarseGridY);
            self.coarseMesh = self.coarseMesh.setBoundaries(2:(2*nX + 2*nY),...
                p_bc_handle, u_bc_handle);
        end
        
        function initialize(self, nData)
            %Initialize model parameters
            %   nFeatures:      number of feature functions
            %   nElements:      number of macro elements
            %   nSCells:        number of cells in S-grid
            
            %Initialize sigma_c
            self.Sigma_c = 1e0*eye(self.gridRF.nCells);
            
            nSX = numel(self.fineGridX); nSY = numel(self.fineGridY);
            if any(self.interpolationMode)
                nSX = nSX + 1; nSY = nSY + 1;
            end
            self.sigma_cf.s0 = ones(nSX*nSY, 1);  %variance field of p_cf
            
            %Initialize hyperparameters
            %this is not working as theta_c is initialized only later
            if strcmp(self.prior_theta_c, 'RVM')
                self.gamma = 1e-4*ones(size(self.theta_c));
            elseif strcmp(self.prior_theta_c, 'VRVM')
                self.gamma = 1e4*ones(size(self.theta_c));
            elseif strcmp(self.prior_theta_c, 'sharedVRVM')
                self.gamma = 1e-4*ones(size(self.theta_c));
            elseif strcmp(self.prior_theta_c, 'adaptiveGaussian')
                self.gamma = 1;
            elseif strcmp(self.prior_theta_c, 'none')
                self.gamma = NaN;
            else
                error('What prior model for theta_c?')
            end
            
            %Initialize parameters of variational approximate distributions
            self.variational_mu{1} = -8.5*ones(1, self.gridRF.nCells);
            self.variational_mu = repmat(self.variational_mu, nData, 1);
            
            self.variational_sigma{1} = .25e-1*ones(1, self.gridRF.nCells);
            self.variational_sigma = repmat(self.variational_sigma, nData, 1);
            
            %Bring variational distribution params in form for
            %unconstrained optimization
            varDistParamsVecInit{1} = [self.variational_mu{1},...
                -2*log(self.variational_sigma{1})];
            self.varDistParamsVec = repmat(varDistParamsVecInit, nData, 1);
        end
        
        function splitRFcells(self, splt_cells)
            %Split cell of random field and initialize corresponding params
            %from coarser random field discretization
            
            nElc = size(self.Sigma_c, 1);
            self.splitted_cells = [self.splitted_cells, splt_cells];
            for splt_cll = splt_cells
                if isnan(self.cell_dictionary(splt_cll))
                    warning('Trying to split an already splitted cell. Skip.')
                else
                    self.gridRF.split_cell(self.gridRF.cells{splt_cll});
                    
                    %extend Sigma_c
                    index = self.cell_dictionary(splt_cll);
                    self.Sigma_c = blkdiag(self.Sigma_c,...
                        self.Sigma_c(index, index)*eye(4));
                    self.Sigma_c(index, :) = [];
                    self.Sigma_c(:, index) = [];
                    
                    %extend theta_c
                    if(~isempty(self.theta_c) && strcmp(self.mode, 'local'))
                        theta_c_mat = reshape(self.theta_c, [], nElc);
                        theta_c = theta_c_mat;
                        theta_c(:, index) = [];
                        theta_c = theta_c(:);
                        self.theta_c=[theta_c;repmat(theta_c_mat(:,index),4,1)];
                    end
                    
                    %extend hyperparameter gamma
                    if(~isempty(self.gamma) && strcmp(self.mode, 'local'))
                        gamma_mat = reshape(self.gamma, [], nElc);
                        gamma = gamma_mat;
                        gamma(:, index) = [];
                        gamma = gamma(:);
                        self.gamma = [gamma; repmat(gamma_mat(:, index), 4, 1)];
                    end
                    
                    %extend Sigma_theta_c
                    if(~isempty(self.Sigma_theta_c)&& strcmp(self.mode,'local'))
                        nFeatures = size(self.Sigma_theta_c, 1)/nElc;
                        Sigma_theta_c_k = self.Sigma_theta_c(...
                            ((index - 1)*nFeatures + 1):(index*nFeatures),...
                            ((index - 1)*nFeatures + 1):(index*nFeatures));
                        self.Sigma_theta_c(((index - 1)*nFeatures + 1):...
                            (index*nFeatures), :) = [];
                        self.Sigma_theta_c(:, ((index - 1)*nFeatures + 1):...
                            (index*nFeatures)) = [];
                        self.Sigma_theta_c = blkdiag(self.Sigma_theta_c, ...
                            Sigma_theta_c_k, Sigma_theta_c_k, ...
                            Sigma_theta_c_k, Sigma_theta_c_k);
                    end
                    
                    %extend variational parameters
                    for n = 1:numel(self.variational_mu)
                        mu_k = self.variational_mu{n}(index);
                        sigma_k = self.variational_sigma{n}(index);
                        self.variational_mu{n}(index) = [];
                        self.variational_sigma{n}(index) = [];
                        self.variational_mu{n} =...
                            [self.variational_mu{n}, mu_k*ones(1, 4)];
                        self.variational_sigma{n} = ...
                            [self.variational_sigma{n}, sigma_k*ones(1,4)];
                        self.varDistParamsVec{n} = [self.variational_mu{n},...
                            -2*log(self.variational_sigma{n})];
                    end
                    
                    %Update cell index dictionary
                    self.cell_dictionary(splt_cll) = nan;
                    self.cell_dictionary((splt_cll + 1):end) = ...
                        self.cell_dictionary((splt_cll + 1):end) - 1;
                    if isnan(self.cell_dictionary(end))
                        self.cell_dictionary = [self.cell_dictionary, ...
                            (self.cell_dictionary(end - 1) + 1):...
                            (self.cell_dictionary(end - 1) + 4)];
                    else
                        self.cell_dictionary = [self.cell_dictionary, ...
                            (self.cell_dictionary(end) + 1):...
                            (self.cell_dictionary(end) + 4)];
                    end
                    %cll_dict = self.cell_dictionary
                end
            end
            self.rf2fem = self.gridRF.map2fine(self.coarseMesh_geometry);
        end
        
        function printCurrentParams(self)
            %Print current model params on screen
            
            if strcmp(self.mode, 'local')
                disp('theta_c: row = feature, column = macro-cell:')
                curr_Sigma_c = full(diag(self.Sigma_c))
                if strcmp(self.prior_theta_c, 'sharedVRVM')
                    curr_gamma = self.gamma(1:...
                        (numel(self.theta_c)/self.gridRF.nCells));
                else
                    curr_gamma = self.gamma;
                end
            else
                curr_Sigma_c = full(diag(self.Sigma_c))
                curr_gamma = self.gamma;
            end
            if(strcmp(self.prior_theta_c, 'sharedVRVM') || ...
                    strcmp(self.prior_theta_c, 'VRVM'))
                activeFeatures =...
                    [find(curr_gamma < 20), curr_gamma(find(curr_gamma < 20))]
            end
            
        end
        
        function fineScaleInterp(self, X)
            %Precompute shape function interp. on every fine scale vertex
            
            nData = numel(X);
            for n = 1:nData
                self.W_cf{n} = shapeInterp(self.coarseMesh, X{n});
            end
        end
        
        function write2file(self, params)
            if ~exist('./data/', 'dir')
                mkdir('./data/');
            end
            
            if strcmp(params, 'elbo')
                filename = './data/elbo';
                elbo = self.elbo;
                save(filename, 'elbo', '-ascii', '-append');
            end
            
            if strcmp(params, 'cell_score')
                filename = './data/cell_score';
                cell_score = self.cell_score';
                save(filename, 'cell_score', '-ascii', '-append');
            end
            
            if strcmp(params, 'cell_score_full')
                filename = './data/cell_score_full';
                cell_score_full = self.cell_score_full';
                save(filename, 'cell_score_full', '-ascii', '-append');
            end
            
            if strcmp(params, 'sigma_cf_score')
                filename = './data/sigma_cf_score';
                sigma_cf_score = self.sigma_cf_score';
                save(filename, 'sigma_cf_score', '-ascii', '-append');
            end
            
            if strcmp(params, 'inv_sigma_cf_score')
                filename = './data/inv_sigma_cf_score';
                inv_sigma_cf_score = self.inv_sigma_cf_score';
                save(filename, 'inv_sigma_cf_score', '-ascii', '-append');
            end
            
            %Optimal params
            %W matrices
            if any(params == 'W')
                for i = 1:numel(self.W_cf)
                    filename = strcat('./data/Wmat', num2str(i));
                    [rowW, colW, valW] = find(self.W_cf{i});
                    WArray = [rowW, colW, valW]';
                    save(filename, 'WArray', '-ascii')
                end
            end
            
            %gamma
            if strcmp(params, 'thetaPriorHyperparam')
                filename = './data/thetaPriorHyperparam';
                thetaPriorHyperparam = self.gamma';
                save(filename, 'thetaPriorHyperparam', '-ascii', '-append');
            end
            
            %theta_c
            if strcmp(params, 'theta_c')
                filename = './data/theta_c';
                tc = self.theta_c';
                save(filename, 'tc', '-ascii', '-append');
            end
            
            %Sigma_theta_c (variance of posterior on theta_c)
            if strcmp(params, 'sigma_theta_c')
                filename = './data/Sigma_theta_c.mat';
                Sigma_theta_c = self.Sigma_theta_c;
                save(filename, 'Sigma_theta_c');
            end
            
            %sigma
            if strcmp(params, 'sigma_c')
                filename = './data/sigma_c';
                sc = full(diag(self.Sigma_c))';
                save(filename, 'sc', '-ascii', '-append');
            end
            
            %S
            if strcmp(params, 'sigma_cf')
                filename = './data/sigma_cf';
                scf = self.sigma_cf.s0';
                onlyFinal = true;
                if onlyFinal
                    save(filename, 'scf', '-ascii');
                else
                    save(filename, 'scf', '-ascii', '-append');
                end
            end
            
            %Parameters of variational distributions on transformed lambda_c
            if strcmp(params, 'vardist')
                varmu = self.variational_mu;
                varsigma = self.variational_sigma;
                save('./data/vardistparams.mat', 'varmu', 'varsigma');
            end
            
        end
        
        function set_summation_matrix(self, X)
            %sets up sum_in_macrocell matrix
            nx = numel(self.fineGridX) + 1;
            ny = numel(self.fineGridY) + 1;
            self.sum_in_macrocell = zeros(self.gridRF.nCells, nx*ny);
            
            kk = 1;
            for k = 1:numel(self.gridRF.cells)
                if isvalid(self.gridRF.cells{k})
                    self.sum_in_macrocell(kk, :) =...
                        self.gridRF.cells{k}.inside(X)';
                    kk = kk + 1;
                end
            end
        end
        
        function compute_elbo(self, N, XMean, XSqMean, X_vtx)
            %General form of elbo allowing model comparison
            %   N:                   number of training samples
            %   XMean, XSqMean:      first and second moments of transformed
            %                        lambda_c
            disp('computing elbo... ')
            assert(~isempty(self.interpolationMode),...
                'Elbo only implemented with fixed dim(U_f)')
            %ONLY VALID IF QoI IS PRESSURE ONLY
            %Short hand notation
            N_dof = numel(self.fineGridX)*numel(self.fineGridY);
            D_c = self.gridRF.nCells;
            aa = self.VRVM_a;
            bb = self.VRVM_b;
            cc = self.VRVM_c;
            dd = self.VRVM_d;
            ee = self.VRVM_e;
            ff = self.VRVM_f;
            D_theta_c = numel(self.theta_c);
            nFeatures = D_theta_c/D_c;
            if strcmp(self.prior_theta_c, 'sharedVRVM')
                D_gamma = nFeatures; %for shared RVM only!
            else
                D_gamma = D_theta_c;
            end
            
            Sigma_lambda_c = XSqMean - XMean.^2;
            %sum over N and macro-cells
            sum_logdet_lambda_c = sum(sum(log(Sigma_lambda_c)));
            
            try
                if self.diag_theta_c
                    logdet_Sigma_theta_ck = zeros(D_c, 1);
                    if(strcmp(self.mode, 'local'))
                        for k = 1:D_c
                            logdet_Sigma_theta_ck(k) = ...
                                sum(log(self.Sigma_theta_c(...
                                ((k-1)*nFeatures + 1):(k*nFeatures))));
                        end
                        logdet_Sigma_theta_c = sum(logdet_Sigma_theta_ck);
                    else
                        logdet_Sigma_theta_c = logdet(self.Sigma_theta_c, 'chol');
                    end
                else
                    logdet_Sigma_theta_ck = zeros(D_c, 1);
                    if(strcmp(self.mode, 'local'))
                        for k = 1:D_c
                            logdet_Sigma_theta_ck(k) = logdet(self.Sigma_theta_c(...
                                ((k - 1)*nFeatures + 1):(k*nFeatures),...
                                ((k - 1)*nFeatures + 1):(k*nFeatures)), 'chol');
                        end
                        logdet_Sigma_theta_c = sum(logdet_Sigma_theta_ck);
                    else
                        logdet_Sigma_theta_c = logdet(self.Sigma_theta_c, 'chol');
                    end
                end
            catch
                logdet_Sigma_theta_c = logdet(self.Sigma_theta_c);
                warning('Sigma_theta_c not pos. def.')
            end
            
            self.elbo = -.5*N*N_dof*log(2*pi) +.5*sum_logdet_lambda_c + ...
                .5*N*D_c + N_dof*(ee*log(ff) + gammaln(self.e) -...
                gammaln(ee)) - self.e*sum(log(self.f)) + D_c*(cc*log(dd) +...
                gammaln(self.c) - gammaln(cc)) -...
                self.c*sum(log(self.d)) + D_gamma*(aa*log(bb) +...
                gammaln(self.a) - gammaln(aa)) - ...
                self.a*sum(log(self.b(1:D_gamma))) + ...
                .5*logdet_Sigma_theta_c + .5*D_theta_c; ...

            self.set_summation_matrix(X_vtx);

            self.cell_score = .5*sum(log(Sigma_lambda_c), 2) - ...
                self.c*log(self.d);
            if strcmp(self.prior_theta_c, 'sharedVRVM')
                self.cell_score = self.cell_score + .5*logdet_Sigma_theta_ck;
            end
            
            f_contribution = - self.e*log(self.f);
            self.cell_score_full = self.cell_score +...
                self.sum_in_macrocell*f_contribution;
            
            self.sigma_cf_score = self.sum_in_macrocell*sqrt(self.sigma_cf.s0);
            self.inv_sigma_cf_score =...
                self.sum_in_macrocell*sqrt(1./self.sigma_cf.s0);
        end
        
        function plot_params(self)
            %Plots the current theta_c
            
            if(self.plt_params && feature('ShowFigureWindows'))
                if ~isfield(self.pParams, 'figParams')
                    self.pParams.figParams = ...
                        figure('units','normalized','outerposition',[0 0 .5 1]);
                end
                
                %short notation
                nSX = numel(self.fineGridX); nSY = numel(self.fineGridY);
                if any(self.interpolationMode)
                    nSX = nSX + 1; nSY = nSY + 1;
                end
                
                colors = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 0 0 0];
                
                sb2 = subplot(3, 2, 2, 'Parent', self.pParams.figParams);
                bar(self.theta_c, 'linewidth', 1, 'Parent', sb2)
                axis(sb2, 'tight');
                axis(sb2, 'fill');
                sb2.XLabel.String = 'component $i$';
                sb2.YLabel.String = '$\theta_{c,i}$';
                
                sb3 = subplot(3, 2, 3, 'Parent', self.pParams.figParams);
                if(~isfield(self.pParams, 'p_sigma') || ...
                        numel(self.pParams.p_sigma) ~= self.gridRF.nCells)
                    %random colors
                    for d = 1:self.gridRF.nCells
                        self.pParams.p_sigma{d} = animatedline('color',...
                            colors(mod(d, 6) + 1, :), 'Parent', sb3);
                    end
                    sb3.XLabel.String = 'iter';
                    sb3.YLabel.String = '$\sigma_k$';
                    sb3.YScale = 'log';
                end
                for d = 1:self.gridRF.nCells
                    addpoints(self.pParams.p_sigma{d}, ...
                        self.EM_iter, full(self.Sigma_c(d, d)));
                end
                axis(sb3, 'tight');
                axis(sb3, 'fill');
                
                sb4 = subplot(3, 2, 4, 'Parent', self.pParams.figParams);
                
                sigma_c_plot = sqrt(self.rf2fem*diag(self.Sigma_c));
                im = imagesc(reshape(sigma_c_plot, self.coarseMesh.nElX,...
                    self.coarseMesh.nElY)', 'Parent', sb4);
                sb4.YDir = 'normal';
                sb4.Title.String = '$\sigma_k$';
                colorbar('Parent', self.pParams.figParams);
                sb4.XTick = [];
                sb4.YTick = [];
                sb4.GridLineStyle = 'none';
                axis(sb4, 'square');
                
                sb5 = subplot(3, 2, 5, 'Parent', self.pParams.figParams);
                if strcmp(self.prior_theta_c, 'sharedVRVM')
                    gam= self.gamma(1:(numel(self.theta_c)/self.gridRF.nCells));
                else
                    gam = self.gamma;
                end
                if ~isfield(self.pParams, 'p_gamma')
                    %random colors
                    for d = 1:numel(gam)
                        self.pParams.p_gamma{d} = animatedline('color',...
                            colors(mod(d, 6) + 1, :), 'Parent', sb5);
                    end
                    sb5.XLabel.String = 'iter';
                    sb5.YLabel.String = '$\gamma$';
                    sb5.YScale = 'log';
                    sb5.YLim(2) = 1e3;
                end
                for d = 1:numel(gam)
                    addpoints(...
                        self.pParams.p_gamma{d}, self.EM_iter, self.gamma(d));
                end
                sp5.YTick = logspace(1e-3, 1e10, 14);
                axis(sb5, 'tight');
                axis(sb5, 'fill');
                
                sb6 = subplot(3, 2, 6, 'Parent', self.pParams.figParams);
                if ~isfield(self.pParams, 'img_S')
                self.pParams.img_S = imagesc(...
                    reshape(sqrt(self.sigma_cf.s0), nSX, nSY), 'Parent', sb6);
                else
                    self.pParams.img_S.CData =...
                        reshape(sqrt(self.sigma_cf.s0),nSX,nSY);
                end
                sb6.Title.String = 'S';
                colorbar('Parent', self.pParams.figParams);
                sb6.GridLineStyle = 'none';
                axis(sb6, 'square');
                sb6.YDir = 'normal';
                sb6.YTick = [];    sb6.XTick = [];
                drawnow;
                
            end
        end
        
        function plotElbo(self, t_tot)
            
            if ~isfield(self.pElbo, 'figure')
                self.pElbo.figure =...
                    figure('units','normalized','outerposition',[0 0 1 1]);
            end
            
            sp = subplot(1, 1, 1, 'Parent', self.pElbo.figure);
            hold(sp, 'on');
            if ~isfield(self.pElbo, 'animatedline')
                self.pElbo.animatedline = animatedline('Parent', sp);
                self.pElbo.animatedline.LineWidth = 2;
                self.pElbo.animatedline.Marker = 'x';
                self.pElbo.animatedline.MarkerSize = 10;
                if nargin > 1
                    sp.XLabel.String = 'comp time';
                else
                    sp.XLabel.String = 'iteration';
                end
                sp.YLabel.String = 'elbo';
            end
            if nargin > 1
                x = t_tot;
            else
                x = self.EM_iter;
            end
            addpoints(self.pElbo.animatedline, x, self.elbo);
            axis(sp, 'tight');
            axis(sp, 'fill');
            drawnow;
        end
        
        function plotCellScores(self)
            %plot adaptive refinement cell scores
            
            if ~isfield(self.pCellScores, 'figure')
                %This is the very first call. Set up things here
                self.pCellScores.figure =...
                    figure('units','normalized','outerposition',[0 0 1 1]);
                
                self.pCellScores.p_active_cells = {};
                self.pCellScores.p_active_cells_S = {};
                self.pCellScores.p_cell_score = {};
            end
                
            %Axes 1
            self.pCellScores.sp1 =...
                subplot(1, 2, 1, 'Parent', self.pCellScores.figure);
            cbp_lambda = colorbar('Parent', self.pCellScores.figure);
            self.pCellScores.sp1.YDir = 'normal';
            axis(self.pCellScores.sp1, 'square');
            
            
            %Axes 2
            self.pCellScores.sp4 =...
                subplot(1, 2, 2, 'Parent', self.pCellScores.figure);
            self.pCellScores.sp4.XLabel.String = 'iteration';
            self.pCellScores.sp4.YLabel.String = 'cell score';
            

            
            imagesc(reshape(self.rf2fem*(-self.cell_score),...
                numel(self.coarseGridX), numel(self.coarseGridY))',...
                'Parent', self.pCellScores.sp1);
            self.pCellScores.sp1.GridLineStyle = 'none';
            self.pCellScores.sp1.XTick = [];
            self.pCellScores.sp1.YTick = [];
            self.pCellScores.sp1.YDir = 'normal';
            self.pCellScores.sp1.Title.String = 'Elbo cell score';
            
            
            
            if(numel(self.cell_score) ~=...
                    numel(self.pCellScores.p_cell_score))
                map = colormap(self.pCellScores.sp4, 'lines');
                for k = (numel(self.pCellScores.p_cell_score) + 1):...
                        numel(self.cell_score)
                    self.pCellScores.p_cell_score{k} = ...
                        animatedline('Parent', self.pCellScores.sp4);
                    self.pCellScores.p_cell_score{k}.LineWidth = .5;
                    self.pCellScores.p_cell_score{k}.Marker = 'x';
                    self.pCellScores.p_cell_score{k}.MarkerSize = 6;
                    self.pCellScores.p_cell_score{k}.Color = map(k, :);
                end
            end
            for k = 1:numel(self.cell_score)
                addpoints(self.pCellScores.p_cell_score{k}, self.EM_iter,...
                    -self.cell_score(k));
            end
            axis(self.pCellScores.sp4, 'tight');
            axis(self.pCellScores.sp4, 'fill');
            
            drawnow;
            
        end
    end
end

