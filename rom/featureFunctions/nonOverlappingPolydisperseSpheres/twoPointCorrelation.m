function [twoPtCorr] = twoPointCorrelation(diskCenters, diskRadii, gridRF,...
    phase, distance)
%Computes the 2-pt. correlation of microstructures with disk exclusions
%pores == voids where fluid can flow
%   diskCenters:         clear
%   diskRadii:           clear
%   gridRF:              random field mesh

diskRadiiSq = diskRadii.^2;
n = 1;
twoPtCorr = zeros(gridRF.nCells, 1);
t_conv = tic;
for cll = gridRF.cells
    if isvalid(cll{1})  %check if cell has been deleted during splitting
        %this loop can be vectorized
        converged = false;
        nSamples = 1;
        while(~converged)
            %x1, x2 are random points in cell
            x1 = (cll{1}.vertices{3}.coordinates -...
                cll{1}.vertices{1}.coordinates).*rand(1, 2) +...
                cll{1}.vertices{1}.coordinates;
            
            %x2 is at distance 'distance' to x1 -- reject if x2 is not in cell
            x2 = [Inf, Inf];
            t_in = tic;
            t_elapsed = 0;
            while(~cll{1}.inside(x2) && t_elapsed < .5)
                phi = 2*pi*rand;
                dx = distance*[cos(phi), sin(phi)];
                x2 = x1 + dx;
                t_elapsed = toc(t_in);
            end
            %Check if the points lie within a circle, i.e. outside domain
            isout1 = any(sum((x1 - diskCenters).^2, 2) < diskRadiiSq');
            isout2 = any(sum((x2 - diskCenters).^2, 2) < diskRadiiSq');
             
            % ~= because we want 'true' for pores and 'false' for solids
            if(isout1 ~= phase && isout2 ~= phase) 
                %both points are in phase 'phase'
                twoPtCorr(n) = (1/nSamples)*((nSamples - 1)*twoPtCorr(n) + 1);
            else
                twoPtCorr(n) = ((nSamples - 1)/nSamples)*twoPtCorr(n);
            end
            mcerr_twoPtCorr = sqrt((twoPtCorr(n) - twoPtCorr(n)^2)/nSamples);
            %Error limits are hard-coded here
            t_elapsed_conv = toc(t_conv);
            if ((mcerr_twoPtCorr/twoPtCorr(n) < 0.001 && nSamples > 5000) ||...
                    t_elapsed_conv > 5)
                converged = true;
            end
            nSamples = nSamples + 1;
        end
        n = n + 1;  %cell index
    end
end

twoPtCorr(twoPtCorr < 0) = 0;

end
