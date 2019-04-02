function [W] = shapeInterp(coarseMesh, x_fine)
%E(e) gives coarse element of fine element e
% [E] = get_coarse_el(domainf.nEl, domainc.nEl, 1:domainf.nEl);

    function [N, E] = shapeFunctionValues2(x)
        %coarse element
        row = sum(coarseMesh.cum_lElY < x(2));
        if row == 0
            row = 1;
        end
        %upper boundary of domain
        col = sum(coarseMesh.cum_lElX < x(1));
        if col == 0
            col = 1;
        end
        %right boundary of domain
        %E is coarse element x is in
        E = (row - 1)*(coarseMesh.nElX) + col;
        
        %shape function values
        N(1)   = (1/coarseMesh.AEl(E))*...
            (x(1) - coarseMesh.lc(E, 2, 1))*(x(2) - coarseMesh.lc(E, 4, 2));
        N(2,1) = -(1/coarseMesh.AEl(E))*...
            (x(1) - coarseMesh.lc(E, 1, 1))*(x(2) - coarseMesh.lc(E, 4, 2));
        N(3)   = (1/coarseMesh.AEl(E))*...
            (x(1) - coarseMesh.lc(E, 1, 1))*(x(2) - coarseMesh.lc(E, 1, 2));
        N(4)   = -(1/coarseMesh.AEl(E))*...
            (x(1) - coarseMesh.lc(E, 2, 1))*(x(2) - coarseMesh.lc(E, 1, 2));
    end

nVertices = size(x_fine, 1);
tic
R = zeros(4*nVertices, 1);
C = zeros(4*nVertices, 1);
Nvec = zeros(4*nVertices, 1);
is = 1;
ie = 4;
%r is the finescale global node number and the row index of W
for r = 1:nVertices
    %coordinate of fine node
    x = [x_fine(r, 1), x_fine(r, 2)];
    [N, E] = shapeFunctionValues2(x);

    %column indices of W, 4 for every local node of E
    c = coarseMesh.globalNodeNumber(E, :);
    R(is:ie) = r;
    C(is:ie) = c;
    Nvec(is:ie) = N;
    is = is + 4;
    ie = ie + 4;
end
W = sparse(R, C, Nvec);
W_assembly_time = toc
end

