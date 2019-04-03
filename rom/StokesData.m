classdef StokesData < handle
    %class for fine scale data of Stokes equation
    
    properties
        %Seldomly changed parameters are to bechanged here
        meshSize = 256
        numberParams = [7.8, 0.2]   %[min, max] pos. number of circ. exclusions
        numberDist = 'logn';
        margins = [0.003, 0.003, 0.003, 0.003]    %[l., u.] margin for imp. phase
        r_params = [-5.23, 0.3]    %[lo., up.] bound on random blob radius
        coordDist = 'GP'
        coordDist_mu = '0.5_0.5'   %only for gauss
        coordDist_cov = 'squaredExponential'
        radiiDist = 'lognGP'
        
        %for GP density distribution
        densityLengthScale = '0.08'
        sigmoidScale = '1.2'
        %for GP on radii
        sigmaGP_r = 0.4
        l_r = 0.05
        origin_rejection = false
        samples
        nSamples
        %base name of file path
        pathname = []
        
        %The properties below are saved as cell arrays (1 cell = 1 sample)
        X       %coordinates of mesh vertices as cell array
        X_interp%coordinates of interpolation mesh, if data is interpolated
        input_bitmap  %bitmap of microstr.; true is pore, false is solid phase
                %onto regular mesh
        P       %pressure at vertices
        U       %velocity at vertices
        bc      %bondary condition coefficients
        cells   %cell-to-vertex map
        cellOfVertex   %Mapping from vertex to cell of S (p_cf variance)
        N_vertices_tot  %total number of vertices in data
        
        %Microstructural data, e.g. centers & radii of circular inclusions
        microstructData
        %Flow boundary conditions; C++ string, only for const. bc here!
        p_bc = '0.0';
        u_bc = {'u_x=1.0-0.0x[1]', 'u_y=1.0-0.0x[0]'}
        %coefficient distribution for randomized bc's
        a_x_m = []      %use fixed bc's if empty
        a_x_s = 1.0
        a_y_m = 0.0
        a_y_s = 1.0
        a_xy_m = 0.0
        a_xy_s = 1.0
        %Design matrix
        designMatrix
        designMatrixSqSum       %precomputation for efficiency
    end
    
    methods
        function self = StokesData(samples)
            %constructor
            self.samples = samples;
            self.nSamples = numel(samples);
        end
        
        function setPathName(self)
            if isempty(self.pathname)
%                 if(strcmp((java.net.InetAddress.getLocalHost.getHostName),...
%                         'workstation1-room0436') || ...
%                         strcmp((java.net.InetAddress.getLocalHost.getHostName),...
%                         'constantin-ThinkPad-T430s'))
%                 self.pathname = strcat(...
%                     '/home/constantin/cluster/python/data/stokesEquation/');
%                 elseif(strcmp((java.net.InetAddress.getLocalHost.getHostName),...
%                         'master'))
%                     %on headnode
%                     self.pathname =...
%                         strcat('/home/constantin/python/data/stokesEquation/');
%                 else
%                     %on subnodes via ethernet
%                     self.pathname = strcat(...
%                         '/home_eth/constantin/python/data/stokesEquation/');
%                 end

                %path for public repository
                self.pathname = strcat('../data/');
                
                if strcmp(self.coordDist, 'GP')
                    self.pathname = char(strcat(self.pathname, 'meshSize=',...
                        num2str(self.meshSize), '/nonOverlappingDisks/margins=',...
                        num2str(self.margins(1)), '_', ...
                        num2str(self.margins(2)), '_', ...
                        num2str(self.margins(3)), '_', ...
                        num2str(self.margins(4)), '/N~', self.numberDist, '/mu=',...
                        sprintf('%.1f', self.numberParams(1)), '/sigma=', ...
                        sprintf('%.1f', self.numberParams(2)), '/x~',...
                        self.coordDist, '/cov=', self.coordDist_cov, '/l=',...
                        self.densityLengthScale, '/sig_scale=',...
                        self.sigmoidScale, '/r~', self.radiiDist));
                elseif strcmp(self.coordDist, 'engineered')
                    self.pathname = char(strcat(self.pathname, 'meshSize=',...
                        num2str(self.meshSize), '/nonOverlappingDisks/margins=',...
                        num2str(self.margins(1)), '_', ...
                        num2str(self.margins(2)), '_', ...
                        num2str(self.margins(3)), '_', ...
                        num2str(self.margins(4)), '/N~', self.numberDist, '/mu=',...
                        sprintf('%.2f', self.numberParams(1)), '/sigma=', ...
                        sprintf('%.1f', self.numberParams(2)), '/x~',...
                        self.coordDist, '/r~', self.radiiDist));
                elseif strcmp(self.coordDist, 'tiles')
                    self.pathname = char(strcat(self.pathname, 'meshSize=',...
                        num2str(self.meshSize), '/nonOverlappingDisks/margins=',...
                        num2str(self.margins(1)), '_', ...
                        num2str(self.margins(2)), '_', ...
                        num2str(self.margins(3)), '_', ...
                        num2str(self.margins(4)), '/N~', self.numberDist, '/mu=',...
                        sprintf('%.1f', self.numberParams(1)), '/sigma=', ...
                        sprintf('%.1f', self.numberParams(2)), '/x~',...
                        self.coordDist, '/r~', self.radiiDist));
                else
                    self.pathname = char(strcat(self.pathname, 'meshSize=',...
                        num2str(self.meshSize), '/nonOverlappingDisks/margins=',...
                        num2str(self.margins(1)), '_', ...
                        num2str(self.margins(2)), '_', ...
                        num2str(self.margins(3)), '_', ...
                        num2str(self.margins(4)), '/N~', self.numberDist, '/mu=',...
                        sprintf('%.2f', self.numberParams(1)), '/sigma=', ...
                        sprintf('%.1f', self.numberParams(2)), '/x~',...
                        self.coordDist, '/mu=', self.coordDist_mu, '/cov=',...
                        self.coordDist_cov, '/r~', self.radiiDist));
                end
                
                
                if strcmp(self.radiiDist, 'lognGP')
                    self.pathname = char(strcat(self.pathname, '/mu=',...
                        sprintf('%.2f', self.r_params(1)), '/sigma=',...
                        num2str(self.r_params(2)), '/sigmaGP_r=',...
                        num2str(self.sigmaGP_r), '/l=', num2str(self.l_r),'/'));
                else
                    self.pathname = char(strcat(self.pathname, '/mu=',...
                        sprintf('%.2f', self.r_params(1)), '/sigma=',...
                        num2str(self.r_params(2)), '/'));
                end
                
                if self.origin_rejection
                    self.pathname=strcat(self.pathname,'origin_rejection=',...
                        num2str(self.origin_rejection), '/');
                end
            end
        end
        
        function [p_bc_handle, u_bc_handle] = returnBCFun(self, n)
            %converts boundary condition strings to handle functions
            %   n:  labels the data point if random boundary conditions are
            %   used
            if nargin < 2
                error('Should only be used with random bc`s. Is data.bc = []?')
            else
                p_bc_handle = str2func(strcat('@(x)', self.p_bc));
                u_bc_handle{1} = @(x) - self.bc{n}(2) - self.bc{n}(3)*x;
                u_bc_handle{2} = @(y) self.bc{n}(1) + self.bc{n}(3)*y;
                u_bc_handle{3} = @(x) self.bc{n}(2) + self.bc{n}(3)*x;
                u_bc_handle{4} = @(y) - self.bc{n}(1) - self.bc{n}(3)*y;
            end
        end
        
        function readData(self, quantities)
            %Reads in Stokes equation data from fenics
            %samples:          samples to load
            %quantities:       identifier for the quantities to load,
            %                  'x' for vertex locations
            %                  'p' for pressure,
            %                  'u' for velocuty,
            %                  'c' for cell-to-vertex map
            %                  'm' for microstructural data
            
            disp('Reading data from disk...')
            self.setPathName;
            
            cellIndex = 1;
            for n = self.samples
                if isempty(self.a_x_m)
                    foldername = char(strcat(self.pathname, 'p_bc=',...
                        self.p_bc, '/', self.u_bc{1}, '_', self.u_bc{2}));
                else
                    foldername = char(strcat(self.pathname, 'p_bc=',...
                        self.p_bc, '/a_x_m=', sprintf('%.1f', self.a_x_m),...
                        '_a_x_s=', sprintf('%.1f', self.a_x_s), 'a_y_m=',...
                        sprintf('%.1f', self.a_y_m), '_a_y_s=',...
                        sprintf('%.1f', self.a_y_s), 'a_xy_m=',...
                        sprintf('%.1f', self.a_xy_m), '_a_xy_s=',...
                        sprintf('%.1f', self.a_xy_s)));
                end
                filename = char(strcat(foldername, '/solution',...
                    num2str(n), '.mat'));
                file = matfile(filename);
                
                if exist(filename, 'file')
                    
                    if contains(quantities, 'x')
                        self.X{cellIndex} = file.x;
                    end
                    
                    if contains(quantities, 'p')
                        self.P{cellIndex} = file.p';
                        if ~isempty(self.a_x_m)
                            self.bc{cellIndex} = file.bc;
                        end
                    end
                    
                    if contains(quantities, 'u')
                        self.U{cellIndex} = file.u;
                    end
                    
                    if contains(quantities, 'c')
                        cellfile = matfile(char(strcat(self.pathname, 'mesh',...
                            num2str(n), '.mat')));
                        self.cells{cellIndex} = cellfile.cells;
                    end
                    
                    if contains(quantities, 'm')
                        datafile = char(strcat(self.pathname,...
                            'microstructureInformation', num2str(n), '.mat'));
                        self.microstructData{cellIndex} = load(datafile);
                    end
                    
                    %preallocation of design matrices
                    cellIndex = cellIndex + 1;
                else
                    self.samples(self.samples == n) = [];
                    self.nSamples = self.nSamples - 1;
                    warning(strcat(filename, 'not found. Skipping sample.'))
                end
            end
            if isempty(self.designMatrix)
                self.designMatrix = cell(cellIndex - 1, 1);
            end
            
            disp('... data loaded to workspace.')
            %This is hard-coded here s.t. it is not forgotten in predictions
            %and computing variance of the data
            % self.removeSpikes('p', 4);
        end
        
        function shiftData(self, interp, quantity, point, value)
            %shifts observable data by constant, s.t. quantity(point) = value
            if nargin < 3
                quantity = 'p';
            end
            if nargin < 4
                point = [0, 0];
            end
            if nargin < 5
                value = 0;
            end
            %closest point to origin is set to p = 0
            if contains(quantity, 'p')
                for n = 1:numel(self.P)
                    if interp
                        dist = sum((self.X_interp{1} - point).^2, 2);
                    else
                        dist = sum((self.X{n} - point).^2, 2);
                    end
                    p_temp = self.P{n};
                    [~, min_dist_i] = min(dist);
                    p_point = p_temp(min_dist_i);
                    p_temp = p_temp - p_point + value;
                    self.P{n} = p_temp;
                end
            else
                error('shifting only implemented for P')
            end
        end
        
        function removeSpikes(self, quantity, spikeLimit)
            %this actually shouldn't be used. It artificially removes spikes 
            %from the fine  scale response
            
            if contains(quantity, 'p')
                %remove spikes from pressure field
                for n = 1:numel(self.P)
                    sorted_p = sort(self.P{n});
                    n_p = numel(sorted_p);
                    cut = round(.05*n_p);
                    mean_p = mean(sorted_p((1 + cut):(end - cut)));
                    std_p = std(sorted_p((1 + cut):(end - cut)));
                    self.P{n}(self.P{n} > mean_p + spikeLimit*std_p) =...
                        mean_p + spikeLimit*std_p;
                    self.P{n}(self.P{n} < mean_p - spikeLimit*std_p) =...
                        mean_p - spikeLimit*std_p;
                end
            end
        end
        
        function interpolate(self, modelParams)
            %Interpolates finescale data onto a regular rectangular grid
            %specified by fineGridX, fineGridY
            
            fineGridX = [0, cumsum(modelParams.fineGridX)];
            fineGridY = [0, cumsum(modelParams.fineGridY)];
            
            if isempty(self.X)
                self.readData('x');
            end
            
            %Specify query grid
            [xq, yq] = meshgrid(fineGridX, fineGridY);
            for n = 1:numel(self.P)
                if ~isempty(self.P)
                    F = scatteredInterpolant(self.X{n}(:, 1),...
                        self.X{n}(:, 2), self.P{n});
                    F.Method = modelParams.interpolationMode;
                    p_interp = F(xq(:), yq(:));
                    
                    %replace original data by interpolated data
                    if ~isempty(modelParams.smoothingParameter)
                        p_interp = reshape(p_interp, numel(fineGridX), ...
                            numel(fineGridY));
                        if modelParams.boundarySmoothingPixels > 0
                            %only smooth boundary
                            p_temp = imgaussfilt(p_interp, ...
                                modelParams.smoothingParameter,...
                                'Padding', 'symmetric');
                            p_interp(1:modelParams.boundarySmoothingPixels,:)...
                                = p_temp(1:...
                                modelParams.boundarySmoothingPixels, :);
                            p_interp((end - ...
                                modelParams.boundarySmoothingPixels):end,:)=...
                                p_temp((end -...
                                modelParams.boundarySmoothingPixels):end,:);
                            p_interp(:,...
                                1:modelParams.boundarySmoothingPixels) = ...
                                p_temp(:,1:modelParams.boundarySmoothingPixels);
                            p_interp(:, (end -...
                                modelParams.boundarySmoothingPixels):end)=...
                                p_temp(:, (end -...
                                modelParams.boundarySmoothingPixels):end);
                        else
                            p_interp = imgaussfilt(p_interp,...
                                modelParams.smoothingParameter,...
                                'Padding', 'symmetric');
                        end
                        p_interp = p_interp(:);
                    end
                    self.P{n} = p_interp;
                end
                
                if ~isempty(self.U)
                    u_interp_x = griddata(self.X{n}(:, 1), self.X{n}(:, 2), ...
                        self.U{n}(1, :), xq(:), yq(:), interpolationMode);
                    u_interp_y = griddata(self.X{n}(:, 1), self.X{n}(:, 2), ...
                        self.U{n}(2, :), xq(:), yq(:), interpolationMode);
                    %replace original data by interpolated data
                    self.U{n} = [];
                    self.U{n} = [u_interp_x, u_interp_y];
                end
            end
            self.X_interp{1} = [xq(:), yq(:)];
        end
        
        function dataVar = computeDataVariance(self, samples, quantity,...
                fineGridX, fineGridY, interpolationMode, smoothingParameter,...
                boundarySmoothingPixels)
            %Computes variance of Stokes data over whole data set
            %keep in mind python indexing for samples
            
            self.samples = samples;
            disp('Reading in data...')
            self.readData(quantity);
            disp('... data read in.')
            
            if isempty(self.X_interp)
                %Data has not yet been interpolated onto a regular grid
                if nargin < 7
                    %no smoothing
                    smoothingParameter = [];
                end
                if nargin < 6
                    interpolationMode = 'cubic';
                end
                if nargin < 5
                    fineGridY = fineGridX;
                end
                self.interpolate(fineGridX, fineGridY, interpolationMode,...
                    smoothingParameter, boundarySmoothingPixels);
            end
            
            %Compute first and second moments
            meanQuantity = 0;
            meanSquaredQuantity = 0;
            for n = 1:numel(samples)
                if strcmp(quantity, 'p')
                    meanQuantity = (1/n)*((n - 1)*meanQuantity + self.P{n});
                    meanSquaredQuantity =...
                        (1/n)*((n - 1)*meanSquaredQuantity + self.P{n}.^2);
                elseif strcmp(quantity, 'u')
                    meanQuantity = (1/n)*((n - 1)*meanQuantity + self.U{n}(:));
                    meanSquaredQuantity =...
                        (1/n)*((n - 1)*meanSquaredQuantity + self.U{n}(:).^2);
                else
                    error('unknown quantity');
                end
            end
            dataVar = meanSquaredQuantity - meanQuantity.^2;
            meanDataVar = mean(dataVar);
        end
        
        function input2bitmap(self, resolution)
            if nargin < 2
                resolution = 256;
            end
            
            if isempty(self.microstructData)
                self.readData('m');
            end
            
            %loop over cirlces
            [xx, yy] = meshgrid(linspace(0, 1, resolution));
            for n = 1:self.nSamples
                r2 = self.microstructData{n}.diskRadii.^2;
                self.input_bitmap{n} = false(resolution);
                
                for nCircle = 1:numel(self.microstructData{n}.diskRadii)
                    self.input_bitmap{n} = self.input_bitmap{n} | ...
                        ((xx - self.microstructData{n}.diskCenters(nCircle, 1)).^2 +...
                        (yy - self.microstructData{n}.diskCenters(nCircle, 2)).^2 <= r2(nCircle));
                end
            end
        end
        
        function countVertices(self)
            self.N_vertices_tot = 0;
            if isempty(self.P)
                self.readData('p');
            end
            for cellIndex = 1:numel(self.P)
                self.N_vertices_tot = self.N_vertices_tot +...
                    numel(self.P{cellIndex});
            end
        end
        
        function vtx2Cell(self, modelParams)
            %Mapping from vertex to cell of rectangular grid specified by
            %grid vectors gridX, gridY

            cumsumX = cumsum(modelParams.fineGridX);
            cumsumX(end) = cumsumX(end) + 1e-12;  %include vertices on boundary
            cumsumY = cumsum(modelParams.fineGridY);
            cumsumY(end) = cumsumY(end) + 1e-12;  %include vertices on boundary
            
            Nx = numel(modelParams.fineGridX);
            
            if(isempty(self.X) && isempty(self.X_interp))
                self = self.readData('x');
            end
            if any(modelParams.interpolationMode)
                if isempty(self.X_interp)
                    self = self.interpolate(modelParams);
                end
                X = self.X_interp;
            else
                X = self.X;
            end
            
            for n = 1:numel(X)
                self.cellOfVertex{n} = zeros(size(X{n}, 1), 1);
                for vtx = 1:size(X{n}, 1)
                    nx = 1;
                    while(X{n}(vtx, 1) > cumsumX(nx))
                        nx = nx + 1;
                    end
                    ny = 1;
                    while(X{n}(vtx, 2) > cumsumY(ny))
                        ny = ny + 1;
                    end
                    self.cellOfVertex{n}(vtx) = nx + (ny - 1)*Nx;
                end
            end
        end
        
        function evaluateFeatures(self, gridRF)
            %Evaluates the feature functions
            if isempty(self.microstructData)
                self.readData('m');
            end
            resolution = 512;
            if isempty(self.input_bitmap)
                self.input2bitmap(resolution);
            end
            
            disp('Evaluating feature functions...');
            
            mData = self.microstructData;
            [indicator, nrow, ncol] = gridRF.indexIndicator(resolution);
            for n = 1:numel(self.samples)
                for k = 1:gridRF.nCells
                    lambdaMat{n, k} = reshape(...
                        self.input_bitmap{n}(indicator{k}),nrow(k), ncol(k));
                end
            end
            dMat = self.designMatrix;
            delta_log = 1;

            grid1x1 = RectangularMesh(1);
            M1 = grid1x1.map2fine(gridRF);
                        
            parPoolInit(self.nSamples);
            parfor n = 1:numel(self.samples)
                %constant 1
                dMat{n} = [dMat{n}, ones(gridRF.nCells, 1)];
                if n == 1
                    dlmwrite('./data/features', 'const', 'delimiter', '');
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% FULL DOMAIN
                %pore fraction
                phi = volumeFractionCircExclusions(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'poreFraction1x1',...
                        'delimiter', '', '-append');
                end
                
                %log pore fraction
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logPoreFraction1x1',...
                        'delimiter', '', '-append');
                end
                
                %sqrt pore fraction
                dMat{n} = [dMat{n}, M1*sqrt(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'sqrtPoreFraction1x1',...
                        'delimiter', '', '-append');
                end
                
                %pore fraction powers (Archie's law, DEM, see Torquato 18.24)
                dMat{n} = [dMat{n}, M1*phi(:).^1.5];
                if n == 1
                    dlmwrite('./data/features', 'PoreFraction^1.5_1x1',...
                        'delimiter', '', '-append');
                end
                
                dMat{n} = [dMat{n}, M1*phi(:).^2];
                if n == 1
                    dlmwrite('./data/features', 'PoreFraction^2_1x1',...
                        'delimiter', '', '-append');
                end
                
                dMat{n} = [dMat{n}, M1*phi(:).^2.5];
                if n == 1
                    dlmwrite('./data/features', 'PoreFraction^2.5_1x1',...
                        'delimiter', '', '-append');
                end
                
                %exp pore fraction
                dMat{n} = [dMat{n}, M1*exp(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'expPoreFraction1x1',...
                        'delimiter', '', '-append');
                end
                
                %log self-consistent approximation (fully ins. spheres)
                dMat{n} = [dMat{n}, M1*log(abs(2*phi(:) - 1) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_SCA1x1',...
                        'delimiter', '', '-append');
                end
                
                %Maxwell approximation
                dMat{n} = [dMat{n}, M1*(phi(:)./(2 - phi(:)))];
                if n == 1
                    dlmwrite('./data/features', 'maxwellApproximation1x1',...
                        'delimiter', '', '-append');
                end
                
                %log Maxwell approximation
                dMat{n} = [dMat{n}, M1*log(phi(:)./(2 - phi(:)) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_maxwellApproximation1x1',...
                        'delimiter', '', '-append');
                end
                
                
                %chord length density
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .005);
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens005_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .0025);
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens0025_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .00125);
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens00125_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .000625);
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens000625_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .0003);
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens0003_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, 0);
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens0_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                
                %interface area
                phi = interfacePerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'interfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log interface area
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logInterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log interface area^1.5
                dMat{n} = [dMat{n}, M1*abs(log(phi(:) + delta_log)).^1.5];
                if n == 1
                    dlmwrite('./data/features', 'abs(logInterfaceArea)^1.5_1x1',...
                        'delimiter', '', '-append');
                end
                
                %square log interface area
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log).^2];
                if n == 1
                    dlmwrite('./data/features', 'log^2InterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %cube log interface area
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log).^3];
                if n == 1
                    dlmwrite('./data/features', 'log^3InterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log^4 interface area
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log).^4];
                if n == 1
                    dlmwrite('./data/features', 'log^4InterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log^1/2 interface area
                dMat{n} = [dMat{n}, M1*abs(log(phi(:) + delta_log)).^.5];
                if n == 1
                    dlmwrite('./data/features', 'log^1/2InterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log^1/3 interface area
                dMat{n} = [dMat{n}, M1*abs(log(phi(:) + delta_log)).^(1/3)];
                if n == 1
                    dlmwrite('./data/features', 'log^1/3InterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log^1/4 interface area
                dMat{n} = [dMat{n}, M1*abs(log(phi(:) + delta_log)).^.25];
                if n == 1
                    dlmwrite('./data/features', 'log^1/4InterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %sqrt interface area
                dMat{n} = [dMat{n}, M1*sqrt(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'sqrtInterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                %interface area ^(1/3)
                dMat{n} = [dMat{n}, M1*phi(:).^(1/3)];
                if n == 1
                    dlmwrite('./data/features', 'InterfaceArea^(1/3)_1x1',...
                        'delimiter', '', '-append');
                end
                
                %interface area ^(1/4)
                dMat{n} = [dMat{n}, M1*phi(:).^(1/4)];
                if n == 1
                    dlmwrite('./data/features', 'InterfaceArea^(1/4)_1x1',...
                        'delimiter', '', '-append');
                end
                
                %interface area ^(1/5)
                dMat{n} = [dMat{n}, M1*phi(:).^(1/5)];
                if n == 1
                    dlmwrite('./data/features', 'InterfaceArea^(1/5)_1x1',...
                        'delimiter', '', '-append');
                end
                
                %square interface area
                dMat{n} = [dMat{n}, M1*phi(:).^2];
                if n == 1
                    dlmwrite('./data/features', 'squareInterfaceArea_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                
                %mean distance between disk edges
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, 'mean',...
                    'edge2edge');
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'meanDist_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log mean distance
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMeanDist_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                %Hagen-Poiseuille?
                %square log mean distance
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log).^2];
                if n == 1
                    dlmwrite('./data/features', 'squareLogMeanDist_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log^3 mean distance
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log).^3];
                if n == 1
                    dlmwrite('./data/features', 'log^3MeanDist_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                %mean distance between disk centers
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, 'mean', 2);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'meanDistCenter_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                %min distance between disk centers
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, 'min', 2);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'minDistCenter_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log min distance
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMinDistCenter_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log min distance squared. Hagen-Poiseuille?
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log).^2];
                if n == 1
                    dlmwrite('./data/features', 'logMinDistCenterSq_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                % lin path
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .25);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath25_1x1',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath25_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .1);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath1_1x1',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath1_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .05);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath05_1x1',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath05_1x1',...
                        'delimiter', '', '-append');
                end
                
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .02);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath02_1x1',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath02_1x1',...
                        'delimiter', '', '-append');
                end
                
                                
                %e_v == porefrac for d == 0
                [~, h_v, poreSizeDens, ~] =...
                    voidNearestSurfaceExclusion(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, 0);

                dMat{n} = [dMat{n}, M1*h_v(:)];
                if n == 1
                    dlmwrite('./data/features', 'h_v0_1x1',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, M1*poreSizeDens(:)];
                if n == 1
                    dlmwrite('./data/features', 'poreSizeProbDens0_1x1',...
                        'delimiter', '', '-append');
                end
                
                
                %log
                dMat{n} =...
                    [dMat{n}, M1*log(h_v(:)+delta_log),...
                    M1*log(poreSizeDens(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_h_v0_1x1',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'log_poreSizeProbDens0_1x1',...
                        'delimiter', '', '-append');
                end
                
                %mean chord length
                phi = meanChordLength(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'meanChordLength_1x1',...
                        'delimiter', '', '-append');
                end
                
                %log mean chord length
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMeanChordLength_1x1',...
                        'delimiter', '', '-append');
                end
                
                %exp mean chord length
                dMat{n} = [dMat{n}, M1*exp(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'expMeanChordLength_1x1',...
                        'delimiter', '', '-append');
                end
                
                %sqrt mean chord length
                dMat{n} = [dMat{n}, M1*sqrt(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'sqrtMeanChordLength_1x1',...
                        'delimiter', '', '-append');
                end
                
                %Radii moments
                phi = momentPerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .2);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', '0.2_moment',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log0.2_moment',...
                        'delimiter', '', '-append');
                end
                
                phi = momentPerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, .5);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', '0.5_moment',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_0.5_moment',...
                        'delimiter', '', '-append');
                end
                
                phi = momentPerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, grid1x1, 1.0);
                dMat{n} = [dMat{n}, M1*phi(:)];
                if n == 1
                    dlmwrite('./data/features', '1.0_moment',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, M1*log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_1.0_moment',...
                        'delimiter', '', '-append');
                end
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %on random field discretization
                %pore fraction
                phi = volumeFractionCircExclusions(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'poreFraction',...
                        'delimiter', '', '-append');
                end
                
                %log pore fraction
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logPoreFraction',...
                        'delimiter', '', '-append');
                end
                
                %sqrt pore fraction
                dMat{n} = [dMat{n}, sqrt(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'sqrtPoreFraction',...
                        'delimiter', '', '-append');
                end
                
                %pore fraction powers (Archie's law, DEM, see Torquato 18.24)
                dMat{n} = [dMat{n}, phi(:).^1.5];
                if n == 1
                    dlmwrite('./data/features', 'PoreFraction^1.5',...
                        'delimiter', '', '-append');
                end
                
                dMat{n} = [dMat{n}, phi(:).^2];
                if n == 1
                    dlmwrite('./data/features', 'PoreFraction^2',...
                        'delimiter', '', '-append');
                end
                
                dMat{n} = [dMat{n}, phi(:).^2.5];
                if n == 1
                    dlmwrite('./data/features', 'PoreFraction^2.5',...
                        'delimiter', '', '-append');
                end
                
                %exp pore fraction
                dMat{n} = [dMat{n}, exp(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'expPoreFraction',...
                        'delimiter', '', '-append');
                end
                
                %log self-consistent approximation (fully ins. spheres)
                dMat{n} = [dMat{n}, log(abs(2*phi(:) - 1) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_SCA',...
                        'delimiter', '', '-append');
                end
                
                %Maxwell approximation
                dMat{n} = [dMat{n}, (phi(:)./(2 - phi(:)))];
                if n == 1
                    dlmwrite('./data/features', 'maxwellApproximation',...
                        'delimiter', '', '-append');
                end
                
                %log Maxwell approximation
                dMat{n} = [dMat{n}, log(phi(:)./(2 - phi(:)) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_maxwellApproximation',...
                        'delimiter', '', '-append');
                end
                
                
                %chord length density
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .005);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens005',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .0025);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens0025',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .00125);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens00125',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .000625);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens000625',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .0003);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens0003',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .00015);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens00015',...
                        'delimiter', '', '-append');
                end
                
                phi = chordLengthDensity(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 0);
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logChordLengthDens0',...
                        'delimiter', '', '-append');
                end
                
                
                
                %interface area
                phi = interfacePerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'interfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log interface area
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logInterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log interface area^1.5
                dMat{n} = [dMat{n}, abs(log(phi(:) + delta_log)).^1.5];
                if n == 1
                    dlmwrite('./data/features', 'abs(logInterfaceArea)^1.5',...
                        'delimiter', '', '-append');
                end
                
                %square log interface area
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^2];
                if n == 1
                    dlmwrite('./data/features', 'log^2InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %cube log interface area
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^3];
                if n == 1
                    dlmwrite('./data/features', 'log^3InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log^4 interface area
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^4];
                if n == 1
                    dlmwrite('./data/features', 'log^4InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log^5 interface area
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^5];
                if n == 1
                    dlmwrite('./data/features', 'log^5InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log^1/2 interface area
                dMat{n} = [dMat{n}, abs(log(phi(:) + delta_log)).^.5];
                if n == 1
                    dlmwrite('./data/features', 'log^1/2InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log^1/3 interface area
                dMat{n} = [dMat{n}, abs(log(phi(:) + delta_log)).^(1/3)];
                if n == 1
                    dlmwrite('./data/features', 'log^1/3InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %log^1/4 interface area
                dMat{n} = [dMat{n}, abs(log(phi(:) + delta_log)).^.25];
                if n == 1
                    dlmwrite('./data/features', 'log^1/4InterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %sqrt interface area
                dMat{n} = [dMat{n}, sqrt(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'sqrtInterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                %interface area ^(1/3)
                dMat{n} = [dMat{n}, phi(:).^(1/3)];
                if n == 1
                    dlmwrite('./data/features', 'InterfaceArea^(1/3)',...
                        'delimiter', '', '-append');
                end
                
                %interface area ^(1/4)
                dMat{n} = [dMat{n}, phi(:).^(1/4)];
                if n == 1
                    dlmwrite('./data/features', 'InterfaceArea^(1/4)',...
                        'delimiter', '', '-append');
                end
                
                %interface area ^(1/5)
                dMat{n} = [dMat{n}, phi(:).^(1/5)];
                if n == 1
                    dlmwrite('./data/features', 'InterfaceArea^(1/5)',...
                        'delimiter', '', '-append');
                end
                
                %square interface area
                dMat{n} = [dMat{n}, phi(:).^2];
                if n == 1
                    dlmwrite('./data/features', 'squareInterfaceArea',...
                        'delimiter', '', '-append');
                end
                
                
                
                %mean distance between disk edges
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'mean',...
                    'edge2edge');
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'meanDist',...
                        'delimiter', '', '-append');
                end
                
                %log mean distance
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMeanDist',...
                        'delimiter', '', '-append');
                end
                
                %max distance between disk edges
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'max',...
                    'edge2edge');
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'maxDist',...
                        'delimiter', '', '-append');
                end
                
                %log max distance
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMaxDist',...
                        'delimiter', '', '-append');
                end
                
                %std distance between disk edges
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'std',...
                    'edge2edge');
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'stdDist',...
                        'delimiter', '', '-append');
                end
                
                %log std distance
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logstdDist',...
                        'delimiter', '', '-append');
                end
                
                %square well potential
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'squareWellPot',...
                    'edge2edge', .01);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'squareWellPot01',...
                        'delimiter', '', '-append');
                end
                
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'squareWellPot',...
                    'edge2edge', .02);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'squareWellPot02',...
                        'delimiter', '', '-append');
                end
                
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'squareWellPot',...
                    'edge2edge', .03);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'squareWellPot03',...
                        'delimiter', '', '-append');
                end
                
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'squareWellPot',...
                    'edge2edge', .04);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'squareWellPot04',...
                        'delimiter', '', '-append');
                end
                
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'squareWellPot',...
                    'edge2edge', .05);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'squareWellPot05',...
                        'delimiter', '', '-append');
                end
                
                
%                 %Lennard-Jones potential
%                 phi = diskDistance(mData{n}.diskCenters,...
%                     mData{n}.diskRadii, gridRF, 'lennardJones',...
%                     'edge2edge', .002);
%                 dMat{n} = [dMat{n}, phi(:)];
%                 if n == 1
%                     dlmwrite('./data/features', 'lennardJones002',...
%                         'delimiter', '', '-append');
%                 end
%                 
%                 phi = diskDistance(mData{n}.diskCenters,...
%                     mData{n}.diskRadii, gridRF, 'lennardJones',...
%                     'edge2edge', .004);
%                 dMat{n} = [dMat{n}, phi(:)];
%                 if n == 1
%                     dlmwrite('./data/features', 'lennardJones004',...
%                         'delimiter', '', '-append');
%                 end
%                 
%                 phi = diskDistance(mData{n}.diskCenters,...
%                     mData{n}.diskRadii, gridRF, 'lennardJones',...
%                     'edge2edge', .006);
%                 dMat{n} = [dMat{n}, phi(:)];
%                 if n == 1
%                     dlmwrite('./data/features', 'lennardJones006',...
%                         'delimiter', '', '-append');
%                 end
%                 
%                 phi = diskDistance(mData{n}.diskCenters,...
%                     mData{n}.diskRadii, gridRF, 'lennardJones',...
%                     'edge2edge', .008);
%                 dMat{n} = [dMat{n}, phi(:)];
%                 if n == 1
%                     dlmwrite('./data/features', 'lennardJones008',...
%                         'delimiter', '', '-append');
%                 end
%                 
%                 phi = diskDistance(mData{n}.diskCenters,...
%                     mData{n}.diskRadii, gridRF, 'lennardJones',...
%                     'edge2edge', .01);
%                 dMat{n} = [dMat{n}, phi(:)];
%                 if n == 1
%                     dlmwrite('./data/features', 'lennardJones01',...
%                         'delimiter', '', '-append');
%                 end
                
                %Hagen-Poiseuille?
                %square log mean distance
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'mean',...
                    'edge2edge');
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^2];
                if n == 1
                    dlmwrite('./data/features', 'squareLogMeanDist',...
                        'delimiter', '', '-append');
                end
                
                %log^3 mean distance
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^3];
                if n == 1
                    dlmwrite('./data/features', 'log^3MeanDist',...
                        'delimiter', '', '-append');
                end
                
                
                %mean distance between disk centers
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'mean', 2);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'meanDistCenter',...
                        'delimiter', '', '-append');
                end
                
                
                %min distance between disk centers
                phi = diskDistance(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 'min', 2);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'minDistCenter',...
                        'delimiter', '', '-append');
                end
                
                %log min distance
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMinDistCenter',...
                        'delimiter', '', '-append');
                end
                
                %log min distance squared. Hagen-Poiseuille?
                dMat{n} = [dMat{n}, log(phi(:) + delta_log).^2];
                if n == 1
                    dlmwrite('./data/features', 'logMinDistCenterSq',...
                        'delimiter', '', '-append');
                end
                
                
                % lin path
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .25);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath25',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath25',...
                        'delimiter', '', '-append');
                end
                
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .1);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath1',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath1',...
                        'delimiter', '', '-append');
                end
                
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .05);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath05',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath05',...
                        'delimiter', '', '-append');
                end
                
                phi = matrixLinealPath(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .02);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPath02',...
                        'delimiter', '', '-append');
                end
                %log
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logLinPath02',...
                        'delimiter', '', '-append');
                end
                
                                
                %e_v == porefrac for d == 0
                [~, h_v, poreSizeDens, ~] =...
                    voidNearestSurfaceExclusion(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 0);

                dMat{n} = [dMat{n}, h_v(:)];
                if n == 1
                    dlmwrite('./data/features', 'h_v0',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, poreSizeDens(:)];
                if n == 1
                    dlmwrite('./data/features', 'poreSizeProbDens0',...
                        'delimiter', '', '-append');
                end
                
                
                %log
                dMat{n} =...
                    [dMat{n}, log(h_v(:)+delta_log),...
                    log(poreSizeDens(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_h_v0',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'log_poreSizeProbDens0',...
                        'delimiter', '', '-append');
                end
                
                %mean chord length
                phi = meanChordLength(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'meanChordLength',...
                        'delimiter', '', '-append');
                end
                
                %log mean chord length
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'logMeanChordLength',...
                        'delimiter', '', '-append');
                end
                
                %exp mean chord length
                dMat{n} = [dMat{n}, exp(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'expMeanChordLength',...
                        'delimiter', '', '-append');
                end
                
                %sqrt mean chord length
                dMat{n} = [dMat{n}, sqrt(phi(:))];
                if n == 1
                    dlmwrite('./data/features', 'sqrtMeanChordLength',...
                        'delimiter', '', '-append');
                end
                
                %Radii moments
                phi = momentPerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .2);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', '0.2_moment',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log0.2_moment',...
                        'delimiter', '', '-append');
                end
                
                phi = momentPerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, .5);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', '0.5_moment',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_0.5_moment',...
                        'delimiter', '', '-append');
                end
                
                phi = momentPerVolume(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, 1.0);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', '1.0_moment',...
                        'delimiter', '', '-append');
                end
                dMat{n} = [dMat{n}, log(phi(:) + delta_log)];
                if n == 1
                    dlmwrite('./data/features', 'log_1.0_moment',...
                        'delimiter', '', '-append');
                end
                
                phi= linPathLengthScale(mData{n}.diskCenters,...
                    mData{n}.diskRadii, gridRF, [0 .02 .05 .1 .25]);
                dMat{n} = [dMat{n}, phi(:)];
                if n == 1
                    dlmwrite('./data/features', 'linPathLengthScale',...
                        'delimiter', '', '-append');
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %features using bitmap input
                phi_tmp = 0;
                for k = 1:gridRF.nCells
                    phi_tmp(k, 1) = distanceTransform(...
                        lambdaMat{n, k}, 'euclidean', false, 'mean');
                    
                    phi_tmp(k, 2) = distanceTransform(...
                        lambdaMat{n, k}, 'euclidean', false, 'var');
                    
                    phi_tmp(k, 3) = distanceTransform(...
                        lambdaMat{n, k}, 'euclidean', false, 'max');
                    
                    phi_tmp(k, 4) = distanceTransform(...
                        lambdaMat{n, k}, 'chessboard', false, 'mean');
                    
                    phi_tmp(k, 5) = distanceTransform(...
                        lambdaMat{n, k}, 'chessboard', false, 'var');
                    
                    phi_tmp(k, 6) = distanceTransform(...
                        lambdaMat{n, k}, 'chessboard', false, 'max');
                    
                    phi_tmp(k, 7) = distanceTransform(...
                        lambdaMat{n, k}, 'cityblock', false, 'mean');
                    
                    phi_tmp(k, 8) = distanceTransform(...
                        lambdaMat{n, k}, 'cityblock', false, 'var');
                    
                    phi_tmp(k, 9) = distanceTransform(...
                        lambdaMat{n, k}, 'cityblock', false, 'max');
                    phi_tmp(k, 10) = gaussLinFilt(lambdaMat{n, k}, nan, 2);
                    phi_tmp(k, 11) = gaussLinFilt(lambdaMat{n, k}, nan, 5);
                    phi_tmp(k, 12) = gaussLinFilt(lambdaMat{n, k}, nan, 10);
                    phi_tmp(k, 13) = isingEnergy(lambdaMat{n, k});
                    phi_tmp(k, 14) = shortestPath(...
                        lambdaMat{n, k}, 'x', resolution, 'quasi-euclidean');
                    phi_tmp(k, 15) = shortestPath(...
                        lambdaMat{n, k}, 'y', resolution, 'quasi-euclidean');
                    phi_tmp(k, 16) = shortestPath(...
                        lambdaMat{n, k}, 'x', resolution, 'cityblock');
                    phi_tmp(k, 17) = shortestPath(...
                        lambdaMat{n, k}, 'y', resolution, 'cityblock');
                end
                dMat{n} = [dMat{n}, phi_tmp];
                if n == 1
                    dlmwrite('./data/features', 'distTransformEuclideanMean',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformEuclideanVar',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformEuclideanMax',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformChessboardMean',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformChessboardVar',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformChessboardMax',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformCityblockMean',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformCityblockVar',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'distTransformCityblockMax',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'gaussLinFilt2',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'gaussLinFilt5',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'gaussLinFilt10',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'isingEnergy',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'shortestPathXquasiEuclidean',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'shortestPathYquasiEuclidean',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'shortestPathXcityblock',...
                        'delimiter', '', '-append');
                    dlmwrite('./data/features', 'shortestPathYcityblock',...
                        'delimiter', '', '-append');
                end
            end
            self.designMatrix = dMat;
            disp('...feature functions evaluated.');
        end
        
        function shapeToLocalDesignMat(self)
            %Sets separate coefficients theta_c for each macro-cell in a single
            %microstructure sample. Don't execute before rescaling/
            %standardization of design Matrix!
            debug = false; %debug mode
            disp(strcat('Using separate feature coefficients theta_c for', ...
                ' each macro-cell in a microstructure...'));
            [nElc, nFeatureFunctions] = size(self.designMatrix{1});
            Phi{1} = zeros(nElc, nElc*nFeatureFunctions);
            nData = numel(self.designMatrix);
            Phi = repmat(Phi, nData, 1);
            
            %Reassemble design matrix
            for n = 1:nData
                for k = 1:nElc
                    Phi{n}(k, ((k - 1)*nFeatureFunctions + 1):...
                        (k*nFeatureFunctions)) = self.designMatrix{n}(k, :);
                end
                Phi{n} = sparse(Phi{n});
            end
            if debug
                firstDesignMatrixBeforeLocal = self.designMatrix{1}
                firstDesignMatrixAfterLocal = full(Phi{1})
                pause
            end
            self.designMatrix = Phi;
            disp('done')
        end
        
        function [featFuncMin, featFuncMax] = computeFeatureFunctionMinMax(self)
            %Computes min/max of feature function outputs over training data,
            %separately for every macro cell
            featFuncMin = self.designMatrix{1};
            featFuncMax = self.designMatrix{1};
            for n = 1:numel(self.designMatrix)
                featFuncMin(featFuncMin > self.designMatrix{n}) =...
                    self.designMatrix{n}(featFuncMin > self.designMatrix{n});
                featFuncMax(featFuncMax < self.designMatrix{n}) =...
                    self.designMatrix{n}(featFuncMax < self.designMatrix{n});
            end
        end
        
        function rescaleDesignMatrix(self, featFuncMin, featFuncMax)
            %Rescale design matrix s.t. outputs are between 0 and 1
            disp('Rescale design matrix...')
            
            if nargin < 3
                [featFuncMin, featFuncMax] = self.computeFeatureFunctionMinMax;
            end
            
            featFuncDiff = featFuncMax - featFuncMin;
            %to avoid irregularities due to rescaling (if every macro cell
            %has the same feature function output). 
            %Like this rescaling maps to 1
            sameOutput = featFuncDiff == 0;
            for n = 1:numel(self.designMatrix)
                self.designMatrix{n} =...
                    (self.designMatrix{n} - featFuncMin)./(featFuncDiff);
                self.designMatrix{n}(sameOutput) = 1;
            end
            
            %Check if design Matrices are real and finite
            self.checkDesignMatrices;
            self.saveNormalization('rescaling', featFuncMin, featFuncMax);
            disp('done')
        end
        
        function checkDesignMatrices(self)
            %Check for finiteness
            for n = 1:numel(self.designMatrix)
                if(~all(all(all(isfinite(self.designMatrix{n})))))
                    dataPoint = n
                    Phi = self.designMatrix{n}
                    pause
                    [coarseElement, featureFunction] = ...
                        ind2sub(size(self.designMatrix{n}),...
                        find(~isfinite(self.designMatrix{n})))
                    warning(strcat('Non-finite design matrix.',...
                        ' Setting non-finite component to 0.'))
                    self.designMatrix{n}(~isfinite(self.designMatrix{n})) = 0;
                elseif(~all(all(all(isreal(self.designMatrix{n})))))
                    warning('Complex feature function output:')
                    dataPoint = n
                    Phi = self.designMatrix{n}
                    pause
                    [coarseElement, featureFunction] =...
                        ind2sub(size(self.designMatrix{n}),...
                        find(imag(self.designMatrix{n})))
                    disp('Ignoring imaginary part...')
                    self.designMatrix{n} = real(self.designMatrix{n});
                end
            end
        end
        
        function saveNormalization(self, type, a, b)
            disp('Saving design matrix normalization...')

            if ~exist('./data')
                mkdir('./data');
            end
            if strcmp(type, 'standardization')
                save('./data/featureFunctionMean', 'a', '-ascii');
                save('./data/featureFunctionSqMean', 'b', '-ascii');
            elseif strcmp(type, 'rescaling')
                save('./data/featureFunctionMin', 'a', '-ascii');
                save('./data/featureFunctionMax', 'b', '-ascii');
            else
                error('Which type of data normalization?')
            end
        end
        
        function [triHandles, meshHandles, pltHandles, figHandle, cb] =...
                plot(self, samples)
            %Plots the fine scale data and returns handles
                        
            %Load data to make sure that not interpolated data is plotted
            self.readData('c');
            self.readData('x');
            self.readData('p');
            self.readData('u');
            self.readData('m');
            
            %figHandle = figure;
            pltIndex = 1;
            resolution = 1024;
            for n = samples
                %figure(figHandle);
                figHandle(n) = figure('units','normalized',...
                    'outerposition',[0 0 1 1]);
                
                %Mesh
                pltHandles(1, pltIndex) = subplot(1, 3, 1);
                
                [xx, yy] = meshgrid(linspace(0, 1, resolution));
                r2 = self.microstructData{n}.diskRadii.^2;
                img = false(resolution);
                
                for k = 1:numel(self.microstructData{n}.diskRadii)
                    img = img | ((xx -...
                        self.microstructData{n}.diskCenters(k, 1)).^2 +...
                        (yy - self.microstructData{n}.diskCenters(k, 2)).^2 ...
                        <= r2(k));
                end
                
                % fig_handle = figure;
                img_handle = imagesc(img, 'Parent', pltHandles(1, pltIndex));
                grid off;
                xticks([0.5 resolution]);
                yticks([0.5 resolution]);
                xticklabels({0, 1});
                yticklabels({0, 1});
                ax = gca;
                ax.YDir = 'normal';
                
                
                
                
                %pressure field
                pltHandles(2, pltIndex) = subplot(1, 3, 2);
                pltHandles(2, pltIndex).Title.String = 'Pressure $P$';
                triHandles(2, pltIndex) =...
                    trisurf(double(self.cells{n}), self.X{n}(:, 1),...
                    self.X{n}(:, 2), self.P{n});
                triHandles(2, pltIndex).LineStyle = 'none';
                axis square;
                axis tight;
                view(3);
                box on;
                grid on;
                xticks([0 .25 .5 .75 1]);
                yticks([0 .25 .5 .75 1]);
                xticklabels({0, '', '', '', 1});
                yticklabels({0, '', '', '', 1});
                pltHandles(2, pltIndex).Title.String = 'Pressure $P$';
%                 cb(1, pltIndex) = colorbar;
%                 cb(1, pltIndex).Label.String = 'pressure p';
%                 cb(1, pltIndex).Label.Interpreter = 'latex';
                
%                 %velocity field (norm)
                u_norm = sqrt(sum(self.U{n}.^2));
%                 pltHandles(2, pltIndex) = subplot(3, N, pltIndex + N);
%                 triHandles(2, pltIndex) = trisurf(self.cells{n},...
%                    self.X{n}(:, 1), self.X{n}(:, 2), u_norm);
%                 triHandles(2, pltIndex).LineStyle = 'none';
%                 axis square;
%                 axis tight;
%                 view(3);
%                 grid off;
%                 box on;
%                 xticks({});
%                 yticks({});
%                 cb(2, pltIndex) = colorbar;
%                 cb(2, pltIndex).Label.String = 'velocity norm $|u|$';
%                 cb(2, pltIndex).Label.Interpreter = 'latex';
                
                %velocity field (norm), 2d
                pltHandles(3, pltIndex) = subplot(1, 3, 3);
                triHandles(3, pltIndex) = trisurf(double(self.cells{n}),...
                   self.X{n}(:, 1), self.X{n}(:, 2), u_norm);
                triHandles(3, pltIndex).LineStyle = 'none';
                axis square;
                axis tight;
                view(2);
                grid off;
                box on;
                xticks({});
                yticks({});
                cb(3, pltIndex) = colorbar;
                pltHandles(3, pltIndex).Title.String = 'velocity norm $|u|$';
                cb(3, pltIndex).Label.Interpreter = 'latex';
                pltIndex = pltIndex + 1;
                
                export_fig(strcat('~/cluster/images/presentationNov18/data_',...
                    num2str(n)), '-png', '-r350');
            end
        end
        
        function [coeff, score, latent, tsquared, explained, mu] =...
                output_pca_analysis(self, fig)
            %performs PCA on the Stokes pressure response
            addpath('./mesh')
            addpath('./FEM')
            
            self.readData('xp');

            temp_params = ModelParams(self.u_bc, self.p_bc);
            %128 should be fine
            temp_params.fineGridX = (1/128)*ones(1, 128);
            temp_params.fineGridY = temp_params.fineGridX;
            
            if isempty(self.X_interp)
                self.interpolate(temp_params);
            end
            temp_P = cell2mat(self.P);
           
            [coeff, score, latent, tsquared, explained, mu] =...
                pca(temp_P');
            
            if nargin < 2
                fig = figure;
            end
            ax = subplot(1,1,1, 'Parent', fig);
            hold on;
            semilogy(cumsum(explained), 'Parent', ax);
            ax.XLabel.String = 'number of components';
            ax.YLabel.String = 'explained variance';
            axis fill
            axis tight
            ax.XLim = [1 30];
        end
    end
end

