function [log_q, d_log_q, Tc] = log_q_n(Xn, Tf_n_minus_mu, W_cf_n, S_cf_n,...
    theta_c, designMatrix,  coarseMesh, transType, transLimits, rf2fem,...
    onlyGrad)

%Xi must be a column vector
if size(Xn, 2) > 1
    Xn = Xn';
end

[lg_p_c, d_lg_p_c] = log_p_c(Xn, designMatrix, theta_c);
[lg_p_cf, d_lg_p_cf, Tc] = log_p_cf(Tf_n_minus_mu, coarseMesh, Xn,...
    W_cf_n, S_cf_n, transType, transLimits, rf2fem, onlyGrad);

log_q = lg_p_cf + lg_p_c;
d_log_q = d_lg_p_c + d_lg_p_cf;

%Finite difference gradient check
FDcheck = false;
if FDcheck
    disp('Gradient check log q_i')
    conductivity =...
        conductivityBackTransform(Xn(1:coarseMesh.nEl), transType, transLimits);
    d = 1e-5;

    latentDim = coarseMesh.nEl;
    gradFD = zeros(latentDim, 1);
    for i = 1:latentDim
        dXi = zeros(latentDim, 1);
        dXi(i) = d;
        conductivityFD = conductivity + conductivity.*dXi(1:coarseMesh.nEl);
        
        [lg_p_c, ~] = log_p_c(Xn + dXi, designMatrix, theta_c);
        [lg_p_cf, ~] = log_p_cf(Tf_n_minus_mu, coarseMesh, Xn + dXi,...
            W_cf_n, S_cf_n, transType, transLimits);
        
        log_qFD = lg_p_cf + lg_p_c;
        gradFD(i) = (log_qFD - log_q)/d;
    end
    
    relgrad = gradFD./d_log_q
    if(any(abs(relgrad - 1) > .1))
        %Note: if d_log_q << d_log_p_c, d_log_p_cf, then this might be due to
        %numerical issues, i.e. FD gradient is unprecise for small log q, it is
        %possible that the FD gradient is unprecise
        conductivity
        conductivityFD
        Xn
        XiFD = Xn + dXi
        d_log_q
        d_lg_p_c
        d_lg_p_cf
        pause 
    end
    
end

end

