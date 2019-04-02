function [log_p, d_log_p, Tc] = log_p_cf(Tf_n_minus_mu, coarseMesh,Xn,...
    W_cf_n, S_cf_n, transType, transLimits, rf2fem, onlyGrad)
%Coarse-to-fine map
%ignore constant prefactor
%log_p = -.5*logdet(S, 'chol') - .5*(Tf - mu)'*(S\(Tf - mu));
%diagonal S

conductivity = conductivityBackTransform(Xn, transType, transLimits);
conductivity = rf2fem*conductivity;

isotropicDiffusivity = true;
if isotropicDiffusivity
    FEMout = heat2d(coarseMesh, conductivity);
else
    D = zeros(2, 2, coarseMesh.nEl);
    for j = 1:coarseMesh.nEl
        D(:, :, j) =  conductivity(j)*eye(2);
    end
    FEMout = heat2d(coarseMesh, D);
end

Tc = FEMout.u;

Tf_n_minus_mu_minus_WTc = Tf_n_minus_mu - W_cf_n*Tc;
%only for diagonal S!
if onlyGrad
    %to avoid unnecessary computation overhead
    log_p = [];
else
    log_p = -.5*(S_cf_n.sumLogS +...
        (S_cf_n.Sinv_vec*(Tf_n_minus_mu_minus_WTc.^2)));
end


if nargout > 1
    %Gradient of FEM equation system w.r.t. conductivities
    
    d_r = FEMgrad(FEMout.naturalTemperatures, coarseMesh);
%     d_rx = d_r;
    if strcmp(transType, 'log')
        %We need gradient of r w.r.t. log conductivities X,
        %multiply each row with resp. conductivity
        %d_rx(1:coarseMesh.nEl, :)= diag(conductivity)*d_r(1:coarseMesh.nEl, :);
        d_rx = d_r.*conductivity;
    elseif strcmp(transType, 'logit')
        %We need gradient w.r.t. x,
        %where x is - log((lambda_up - lambda_lo)/(lambda - lambda_lo) - 1)
        X = conductivityTransform(conductivity, transType, transLimits);
        dLambda_dX = (transLimits(2) - transLimits(1))./(exp(X) + 2 + exp(-X));
        d_rx(1:coarseMesh.nEl, :) = diag(dLambda_dX)*d_r(1:coarseMesh.nEl, :);
    elseif strcmp(transType, 'log_lower_bound')
        %transformation is X = log(Lambda - lambda_lo)
        dLambda_dX = conductivity - transLimits(1);
        d_rx(1:coarseMesh.nEl, :) = diag(dLambda_dX)*d_r(1:coarseMesh.nEl, :);
    elseif strcmp(transType, 'square')
        %We need gradient of r w.r.t. sqrt conductivities X. d/dX = 2X d/dlambda
        d_rx(1:coarseMesh.nEl, :) = diag(2*Xn(1:coarseMesh.nEl))*...
            d_r(1:coarseMesh.nEl, :);
    else
        error('Unknown conductivity transformation')
    end
    adjoints = get_adjoints(FEMout.globalStiffness,...
        S_cf_n.WTSinv, coarseMesh, Tf_n_minus_mu_minus_WTc);
    d_log_p = - d_rx*adjoints;
    
    %We need the gradient w.r.t. degrees of freedom of discretized rand. field
    d_log_p = rf2fem'*d_log_p;

    
    %Finite difference gradient check
    FDcheck = false;
    if FDcheck
        disp('Gradient check log p_cf')
        if isempty(log_p)
            log_p = -.5*(S_cf_n.sumLogS +...
                (S_cf_n.Sinv_vec*(Tf_n_minus_mu_minus_WTc.^2)));
        end
        d = 1e-8;
        FDgrad = zeros(coarseMesh.nEl, 1);
        for e = 1:coarseMesh.nEl
            conductivityFD = conductivity;
            conductivityFD(e) = conductivityFD(e) + d;
            
            if isotropicDiffusivity
                FEMoutFD = heat2d(coarseMesh, conductivityFD);
            else
                DFD = zeros(2, 2, coarseMesh.nEl);
                for j = 1:coarseMesh.nEl
                    DFD(:, :, j) =  conductivityFD(j)*eye(2);
                end
                FEMoutFD = heat2d(coarseMesh, DFD);
            end

            TcFD = FEMoutFD.Tff';
            TcFD = TcFD(:);
            
            WTcFD = W_cf_n*TcFD;
            log_pFD = -.5*(S_cf_n.sumLogS + (Tf_n_minus_mu - WTcFD)'*...
                (S_cf_n.Sinv_vec.*(Tf_n_minus_mu - WTcFD)));
            if strcmp(transType, 'log')
                FDgrad(e) = conductivity(e)*(log_pFD - log_p)/d;
            elseif strcmp(transType, 'logit')
                FDgrad(e) = dLambda_dX(e)*(log_pFD - log_p)/d;
            elseif strcmp(transType, 'log_lower_bound')
                FDgrad(e) = dLambda_dX(e)*(log_pFD - log_p)/d;
            else
                error('Unknown conductivity transformation')
            end
        end
        
        FDgrad = rf2fem'*FDgrad;
        relgrad = FDgrad./d_log_p
%         plot(1:numel(FDgrad), FDgrad, 1:numel(FDgrad), d_log_p)
%         axis square
%         drawnow
        pause(.1)
        if(norm(relgrad - 1) > 1e-1)
            log_p
            log_pFD
            d_log_p
            FDgrad
            diff = log_pFD - log_p
        end
    end %FD gradient check
    
    
end

end

