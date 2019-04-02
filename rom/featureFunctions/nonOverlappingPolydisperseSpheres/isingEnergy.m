function [E] = isingEnergy(lambdakMat)
%Compute the ising energy of the conductivity field lambda as a feature function

[nr, nc] = size(lambdakMat);

%not efficient but no bottleneck
E = 0;
for c = 1:nc
    for r = 1:nr
        if(r > 1)
            if(lambdakMat(r, c) == lambdakMat(r - 1, c))
                E = E + 1;
            else
                E = E - 1;
            end
        else
            %r == 1; periodic boundary conditions
            if(lambdakMat(1, c) == lambdakMat(nr, c))
                E = E + 1;
            else
                E = E - 1;
            end
        end
        if(r < nr)
            if(lambdakMat(r, c) == lambdakMat(r + 1, c))
                E = E + 1;
            else
                E = E - 1;
            end
        else
            %r == nr; periodic boundary conditions
            if(lambdakMat(nr, c) == lambdakMat(1, c))
                E = E + 1;
            else
                E = E - 1;
            end
        end
        if(c > 1)
            if(lambdakMat(r, c) == lambdakMat(r, c - 1))
                E = E + 1;
            else
                E = E - 1;
            end
        else
            %c == 1; periodic boundary conditions
            if(lambdakMat(r, 1) == lambdakMat(r, nc))
                E = E + 1;
            else
                E = E - 1;
            end
        end
        if(c < nc)
            if(lambdakMat(r, c) == lambdakMat(r, c + 1))
                E = E + 1;
            else
                E = E - 1;
            end
        else
            %c == nc; periodic boundary conditions
            if(lambdakMat(r, nc) == lambdakMat(r, 1))
                E = E + 1;
            else
                E = E - 1;
            end
        end
    end
end

end

