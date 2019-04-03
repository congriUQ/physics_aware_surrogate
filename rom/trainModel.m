%% This is the main training file for the Darcy-type ROM for Stokes equation
%% Preamble:

rehash
restoredefaultpath
clear
addpath('./featureFunctions/nonOverlappingPolydisperseSpheres')
addpath('./mesh')
addpath('./aux')
addpath('./comp')
addpath('./FEM')
addpath('rom')
addpath('./VI')
%No badly conditioned warnings -- turn on if desired
warning('off', 'MATLAB:nearlySingularMatrix');

startup;    %for proper graphics configuration
%random number seed based on time
rng('shuffle');


%% Initialization
%Which data samples for training?
nTrain = 16;
nStart = 0;
samples = nStart:(nTrain - 1 + nStart);
loadParams = false;     %load parameters from previous run?

rom = StokesROM;

rom.trainingData = StokesData(samples);
rom.trainingData.readData('px');
rom.trainingData.countVertices();

if loadParams
    disp('Loading modelParams...')
    load('./data/modelParams.mat');
    rom.modelParams = copy(modelParams);
    if any(rom.modelParams.interpolationMode)
        rom.trainingData.interpolate(rom.modelParams);
        rom.modelParams.fineScaleInterp(rom.trainingData.X_interp);
        interp = true;
    else
        rom.modelParams.fineScaleInterp(rom.trainingData.X);%for W_cf
        interp = false;
    end
    disp('... modelParams loaded.')
else
    rom.modelParams = ModelParams(rom.trainingData.u_bc, rom.trainingData.p_bc);
    rom.modelParams.initialize(nTrain);
    if any(rom.modelParams.interpolationMode)
        rom.trainingData.interpolate(rom.modelParams);
        rom.modelParams.fineScaleInterp(rom.trainingData.X_interp);
        interp = true;
    else
        rom.modelParams.fineScaleInterp(rom.trainingData.X);%for W_cf
        interp = false;
    end
end

%do not remove! if no cell is splitted, pass empty array
rom.modelParams.splitRFcells([]);

%Parameters from previous runs are deleted here
if exist('./data/', 'dir')
    rmdir('./data', 's'); %delete old params
end
mkdir('./data/');

rom.trainingData.shiftData(interp, 'p'); %shifts p to 0 at origin
rom.trainingData.vtx2Cell(rom.modelParams);

sw0_mu = 6e-4;
sw0_sigma = 6e-5;
sw_decay = .995; %decay factor per iteration
VI_t = [60*ones(1, 1), 20*ones(1, 10), 10];

split_schedule = [];
if isempty(split_schedule)
    nSplits = 0;
else
    nSplits = numel(split_schedule);
end
tic_tot = tic;
for split_iter = 1:(nSplits + 1)
    
    clear XMean XSqMean;
    %delete design matrices so that they can be recomputed
    rom.trainingData.designMatrix = cell(1, nTrain);
    rom.trainingData.evaluateFeatures(rom.modelParams.gridRF);
    
    if strcmp(rom.modelParams.normalization, 'rescale')
        rom.trainingData.rescaleDesignMatrix;
    end
    
    if strcmp(rom.modelParams.mode, 'local')
        rom.trainingData.shapeToLocalDesignMat;
    end
        
    %theta_c must be initialized after design matrices exist
    if isempty(rom.modelParams.theta_c)
        disp('Initializing theta_c...')
        rom.modelParams.theta_c =...
            0*ones(size(rom.trainingData.designMatrix{1}, 2), 1);
        disp('...done.')
    end
    
    %Step width for stochastic optimization in VI
    nRFc = rom.modelParams.gridRF.nCells;
    sw = (sw_decay^rom.modelParams.EM_iter)*...
        [sw0_mu*ones(1, nRFc), sw0_sigma*ones(1, nRFc)];
    sw_min = 3e-1*[sw0_mu*ones(1, nRFc), sw0_sigma*ones(1, nRFc)];
    sw(sw < sw_min) = sw_min(sw < sw_min);
        
    if isempty(rom.trainingData.a_x_m)
        %fixed bc's
        coarseMesh = rom.modelParams.coarseMesh;
        coarseMesh = coarseMesh.shrink;
    else
        nX = length(rom.modelParams.coarseGridX);
        nY = length(rom.modelParams.coarseGridY);
        for n = 1:nTrain
            %random bc's
            coarseMesh{n} = rom.modelParams.coarseMesh;
            [p_bc_handle, u_bc_handle] = rom.trainingData.returnBCFun(n); 
            coarseMesh{n} = coarseMesh{n}.setBoundaries(2:(2*nX + 2*nY),...
                p_bc_handle, u_bc_handle);
            coarseMesh{n} = coarseMesh{n}.shrink;
        end
    end
    
    
    %% Actual training phase:
    converged = false;
    ppool = parPoolInit(nTrain);
    pend = 0;
    rom.modelParams.EM_iter_split = 0;
    %reset epoch, important for splitting
    rom.modelParams.epoch = 0;
    while ~converged
        
        %% Setting up a handle to the distribution q_n - this transfers less 
        %data in parfor loops
        for n = 1:nTrain
            P_n_minus_mu = rom.trainingData.P{n};
            if any(rom.modelParams.interpolationMode)
                W_cf_n = rom.modelParams.W_cf{1};
                S_n = rom.modelParams.sigma_cf.s0;
            else
                W_cf_n = rom.modelParams.W_cf{n};
                S_n = rom.modelParams.sigma_cf.s0(...
                    rom.trainingData.cellOfVertex{n});
            end
            %S_n is a vector of variances at vertices
            S_cf_n.sumLogS = sum(log(S_n));
            S_cf_n.Sinv_vec = (1./S_n)';%row vector
            Sinv = sparse(1:length(S_n), 1:length(S_n), S_cf_n.Sinv_vec);
            S_cf_n.WTSinv = (Sinv*W_cf_n)';
            
            tc.theta = rom.modelParams.theta_c;
            tc.Sigma = rom.modelParams.Sigma_c;
            tc.SigmaInv = inv(tc.Sigma);
            
            Phi_n = rom.trainingData.designMatrix{n};
            
            rf2fem = rom.modelParams.rf2fem;
            transType = rom.modelParams.diffTransform;
            transLimits = rom.modelParams.diffLimits;
            if isempty(rom.trainingData.a_x_m)
                cm = coarseMesh;
            else
                cm = coarseMesh{n};
            end
            lg_q{n} = @(Xi) log_q_n(Xi, P_n_minus_mu, W_cf_n, S_cf_n, tc,...
                Phi_n, cm, transType, transLimits, rf2fem, true);
        end
        varDistParamsVec = rom.modelParams.varDistParamsVec;
        
        if(rom.modelParams.epoch > 0 && ~loadParams)
            %Sequentially update N_threads qi's at a time, then perform M-step
            pstart = pend + 1;
            if pstart > nTrain
                pstart = 1;
                rom.modelParams.epoch = rom.modelParams.epoch + 1;
                %Gradually reduce VI step width
                sw = sw_decay*sw;
                sw(sw < sw_min) = sw_min(sw < sw_min);
            end
            pend = pstart + ppool.NumWorkers - 1;
            if pend > nTrain
                pend = nTrain;
            elseif pend < pstart
                pend = pstart;
            end
        else
            pstart = 1;
            pend = nTrain;
            rom.modelParams.epoch = rom.modelParams.epoch + 1;
            loadParams = false; %first iteration after loading of params needs
                                %to cycle through all data
        end
        
        disp('Variational Inference...')
        if(rom.modelParams.epoch + 1 <= numel(VI_t))
            t = VI_t(rom.modelParams.epoch + 1);
        else
            t = VI_t(end);
        end
        ticBytes(gcp);
        tic
        parfor n = pstart:pend
            %Finding variational approximation to q_n
            [varDistParams{n}, varDistParamsVec{n}] = efficientStochOpt(...
                varDistParamsVec{n}, lg_q{n}, 'diagonalGauss', sw, nRFc, t);
        end
        tocBytes(gcp)
        VI_time = toc
        rom.modelParams.varDistParamsVec = varDistParamsVec;
        disp('... VI done.')
        
        for n = pstart:pend
            %Compute expected values under variational approximation
            rom.modelParams.variational_mu{n} = varDistParams{n}.mu;
            rom.modelParams.variational_sigma{n} = varDistParams{n}.sigma;
            XMean(:, n) = varDistParams{n}.mu';
            XSqMean(:, n) = varDistParams{n}.XSqMean;
            
            P_n_minus_mu = rom.trainingData.P{n};
            if any(rom.modelParams.interpolationMode)
                W_cf_n = rom.modelParams.W_cf{1};
            else
                W_cf_n = rom.modelParams.W_cf{n};
            end
            if isempty(rom.trainingData.a_x_m)
                cm = coarseMesh;
            else
                cm = coarseMesh{n};
            end
            p_cf_expHandle_n = @(X) sqMisfit(X, transType, transLimits,...
                cm, P_n_minus_mu, W_cf_n, rf2fem);
            %Expectations under variational distributions
            p_cf_exp =...
                mcInference(p_cf_expHandle_n,'diagonalGauss', varDistParams{n});
            sqDist{n} = p_cf_exp;
        end
        
        
        %M-step: determine optimal parameters given the sample set
        tic
        disp('M-step...')
        rom.M_step(XMean, XSqMean, sqDist)
        disp('...M-step done.')
        M_step_time = toc

        rom.modelParams.compute_elbo(nTrain, XMean, XSqMean,...
            rom.trainingData.X_interp{1});
        elbo = rom.modelParams.elbo
        
        if(~mod(rom.modelParams.EM_iter_split - 1, 3) && nSplits)
            rom.modelParams.active_cells = rom.findMeshRefinement(false)';
            activeCells = rom.modelParams.active_cells
            filename = './data/activeCells';
            save(filename, 'activeCells', '-ascii', '-append');
        end
        
        %Print some output
        Lambda_eff1_mode = conductivityBackTransform(...
            rom.trainingData.designMatrix{1}*rom.modelParams.theta_c,...
            transType, transLimits)
        rom.modelParams.printCurrentParams;
        
        plt = true;
        if(plt && feature('ShowFigureWindows'))
            disp('Plotting...')
            t_plt = tic;
            %plot parameters
            rom.modelParams.plot_params();
            %plot modal lambda_c and corresponding -training- data reconstruction
            rom.plotCurrentState(0, transType, transLimits);
            %plot elbo vs. training iteration
            t_tot = toc(tic_tot)
            rom.modelParams.plotElbo(t_tot);
            %Plot adaptive refinement cell scores
            rom.modelParams.plotCellScores();
            disp('...plotting done. Plotting time:')
            t_plt = toc(t_plt)
        end
        
        
        %write parameters to disk to investigate convergence
        rom.modelParams.write2file('thetaPriorHyperparam');
        rom.modelParams.write2file('theta_c');
        rom.modelParams.write2file('sigma_c');
        rom.modelParams.write2file('elbo');
        rom.modelParams.write2file('cell_score');
        rom.modelParams.write2file('cell_score_full');
        rom.modelParams.write2file('sigma_cf_score');
        rom.modelParams.write2file('inv_sigma_cf_score');
        if rom.modelParams.epoch > rom.modelParams.max_EM_epochs(...
                min([split_iter, numel(rom.modelParams.max_EM_epochs)]))
            converged = true;
        end
        if(~mod(rom.modelParams.EM_iter, 1) || converged)
            tic
            modelParams = copy(rom.modelParams);
            modelParams.pParams = [];
            modelParams.pElbo = [];
            modelParams.pCellScores = [];
            %Save modelParams after every iteration
            disp('Saving modelParams...')
            save('./data/modelParams.mat', 'modelParams', '-v7.3');
            disp('...modelParams saved.')
            save_time = toc
        end
        save('./data/XMean', 'XMean');
        save('./data/XSqMean', 'XSqMean');
        XMeanVec = XMean(:)';
        save('./data/XMeanVec', 'XMeanVec', '-ascii', '-append');
        XSqMeanVec = XSqMean(:)';
        save('./data/XSqMeanVec', 'XSqMeanVec', '-ascii', '-append');
        epoch = rom.modelParams.epoch
    end
    
    if split_iter < (nSplits + 1)
        disp('splitting cell...')
        refinement_objective = 'active_cells';
        if strcmp(refinement_objective, 'active_cells_S')
            rom.modelParams.active_cells_S = rom.findMeshRefinement(true)';
            activeCells_S = rom.modelParams.active_cells_S;
            filename = './data/activeCells_S';
            save(filename, 'activeCells_S', '-ascii', '-append');
            [~, cell_indices_pde] =...
                sort(rom.modelParams.active_cells_S, 'descend');
        elseif strcmp(refinement_objective, 'active_cells')
            rom.modelParams.active_cells = rom.findMeshRefinement(false)';
            activeCells = rom.modelParams.active_cells;
            filename = './data/activeCells';
            save(filename, 'activeCells', '-ascii', '-append');
            [~, cell_indices_pde] =...
                sort(rom.modelParams.active_cells, 'descend');
        elseif strcmp(refinement_objective, 'full_elbo_score')
            [~, cell_indices_pde] = sort(rom.modelParams.cell_score_full);
        elseif strcmp(refinement_objective, 'sigma_cf')
            [~, cell_indices_pde] =...
                sort(rom.modelParams.sigma_cf_score, 'descend');
        elseif strcmp(refinement_objective, 'inv_sigma_cf')
            [~, cell_indices_pde] =...
                sort(rom.modelParams.inv_sigma_cf_score, 'descend');
        elseif strcmp(refinement_objective, 'reduced_elbo_score')
            [~, cell_indices_pde] = sort(rom.modelParams.cell_score);
        elseif strcmp(refinement_objective, 'random')
            cell_indices_pde = randperm(numel(rom.modelParams.cell_score));
        end
        
        splitable = false;
        split_attempt = 1;
        if isempty(split_schedule)
            while ~splitable
                cell_index = find(rom.modelParams.cell_dictionary ==...
                    cell_indices_pde(split_attempt))
                if rom.modelParams.coarseMesh.AEl(...
                        cell_indices_pde(split_attempt)) <=...
                        .25*rom.modelParams.gridRF.cells{cell_index}.surface
                    rom.modelParams.splitRFcells([cell_index]);
                    splitable = true;
                    disp('...cell splitted.')
                else
                    warning('Cell not splitable. Trying to split second choice...')
                end
                split_attempt = split_attempt + 1;
            end
        else
            rom.modelParams.splitRFcells([split_schedule(split_iter)]);
            disp('...cell splitted.')
        end
        %needs to be recomputed
        rom.trainingData.designMatrixSqSum = [];
    end
end