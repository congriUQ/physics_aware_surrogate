function [log_p, d_log_p, data] = log_p_c(Xq, Phi, theta_c)
%Probabilistic mapping from fine to coarse heat conductivity
%   Xq:         Effective log conductivity vector
%   Phi:        Design Matrix
%   theta_c:    distribution parameters

mu  = Phi*theta_c.theta;    %mean

%ignore constant prefactor
% log_p = - size(Xq, 1)*log(sigma) - (1/(2*sigma^2))*(Xq - mu)'*(Xq - mu);
%Diagonal covariance matrix Sigma
% log_p = - .5*sum(log(diag(theta_c.Sigma)))
%         - .5*(Xq - mu)'*theta_c.SigmaInv*(Xq - mu);
log_p = - .5*logdet(theta_c.Sigma) - .5*(Xq - mu)'*theta_c.SigmaInv*(Xq - mu);

if nargout > 1
    d_log_p = theta_c.SigmaInv*(mu - Xq);
    
    %Finite difference gradient check
    FDcheck = false;
    if FDcheck
        disp('Gradient check log p_c')
        d = 1e-5;
        d_log_pFD = 0*Xq;
        for i = 1:size(Xq, 1)
            dXq = 0*Xq;
            dXq(i) = d;
            d_log_pFD(i) = (-.5*(Xq + dXq - mu)'*...
                (theta_c.SigmaInv*(Xq + dXq - mu)) - log_p)/d;
        end 
%         d_log_pFD
%         d_log_p
        relGrad = d_log_pFD./d_log_p
    end
end

%dummy
if nargout > 2
    data = 0;
end
    
end

