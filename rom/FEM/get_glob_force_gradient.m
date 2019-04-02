function [F] = get_glob_force_gradient(domain, d_k, el)
%Assemble global force vector

%nested function for performance
    function f2 = get_loc_force_gradient2(domain, d_k)
        f2 = zeros(4, 1);
        %Boundary value temperature of element e
        Tb = zeros(4, 1);
        Tbflag = false;
        
        globNode = domain.globalNodeNumber(el, :);
        for i = 1:4
            if(~isnan(domain.essentialTemperatures(globNode(i))))
                Tb(i) = domain.essentialTemperatures(globNode(i));
                Tbflag = true;
            end
        end
        if(Tbflag)
            f2 = -d_k*Tb;
        end
        
    end


F = zeros(domain.nEq,1);
f = get_loc_force_gradient2(domain, d_k);

ln = 1:4;
eqnTemp = domain.lm(el, :);
eqn = eqnTemp(eqnTemp > 0);
ln = ln(eqnTemp > 0);
F(eqn) = F(eqn) + f(ln);

end