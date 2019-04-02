classdef RectangularMesh < Mesh
    %Class for mesh consisting of generic rectangles

    
    methods
        function self = RectangularMesh(gridX, gridY)
            %constructor
            %gridX, gridY are vectors of edge lengths
            
            %Generate mesh according to grid vectors
            if nargin > 0
                if nargin == 1
                    gridY = gridX;
                end
                
                %Create vertices
                x_coord = cumsum([0, gridX]);
                y_coord = cumsum([0, gridY]);
                for y = y_coord
                    for x = x_coord
                        self.create_vertex([x, y]);
                    end
                end
                
                %Create edges
                nx = numel(gridX) + 1;
                ny = numel(gridY) + 1;
                for y = 0:(ny - 1)
                    for x = 1:nx
                        if self.vertices{x + y*nx}.coordinates(1) > 0
                            self.create_edge(self.vertices{x + y*nx - 1}, ...
                                self.vertices{x + y*nx});
                        end
                    end
                end
                for y = 0:(ny - 1)
                    for x = 1:nx
                        if self.vertices{x + y*nx}.coordinates(2) > 0
                            self.create_edge(self.vertices{x + (y - 1)*nx}, ...
                                self.vertices{x + y*nx});
                        end
                    end
                end
                
                %Create cells
                nx = nx - 1;
                ny = ny - 1;
                n = 1;  %cell index
                for y = 1:ny
                    for x = 1:nx
                        vtx = {self.vertices{x + (y - 1)*(nx + 1)}, ...
                            self.vertices{x + (y - 1)*(nx + 1) + 1}, ...
                            self.vertices{x + y*(nx + 1) + 1}, ...
                            self.vertices{x + y*(nx + 1)}};
                        edg= {self.edges{n}, self.edges{nx*(ny + 1) + n + y},...
                            self.edges{n + nx},...
                            self.edges{nx*(ny + 1) + n + y - 1}};
                        self.create_cell(vtx, edg);
                        n = n + 1;
                    end
                end
            end
        end
        
        function self = split_cell(self, cll)
            %Splits the rectangular cell cll

            %Create new vertices
            new_vertices{1} = self.create_vertex(cll.centroid);
            lo_coord = .5*(cll.edges{1}.vertices{1}.coordinates + ...
                cll.edges{1}.vertices{2}.coordinates);
            new_vertices{2} = self.create_vertex(lo_coord);
            r_coord = .5*(cll.edges{2}.vertices{1}.coordinates + ...
                cll.edges{2}.vertices{2}.coordinates);
            new_vertices{3} = self.create_vertex(r_coord);
            up_coord = .5*(cll.edges{3}.vertices{1}.coordinates + ...
                cll.edges{3}.vertices{2}.coordinates);
            new_vertices{4} = self.create_vertex(up_coord);
            le_coord = .5*(cll.edges{4}.vertices{1}.coordinates + ...
                cll.edges{4}.vertices{2}.coordinates);
            new_vertices{5} = self.create_vertex(le_coord);
            
            %Create new edges. Go around old cell, then middle cross
            %go around
            new_edges{1} = self.create_edge(cll.vertices{1}, new_vertices{2});
            new_edges{2} = self.create_edge(new_vertices{2}, cll.vertices{2});
            new_edges{3} = self.create_edge(cll.vertices{2}, new_vertices{3});
            new_edges{4} = self.create_edge(new_vertices{3}, cll.vertices{3});
            new_edges{5} = self.create_edge(cll.vertices{3}, new_vertices{4});
            new_edges{6} = self.create_edge(new_vertices{4}, cll.vertices{4});
            new_edges{7} = self.create_edge(cll.vertices{4}, new_vertices{5});
            new_edges{8} = self.create_edge(new_vertices{5}, cll.vertices{1});
            %middle cross
            new_edges{9} = self.create_edge(new_vertices{1}, new_vertices{2});
            new_edges{10} = self.create_edge(new_vertices{1}, new_vertices{3});
            new_edges{11} = self.create_edge(new_vertices{1}, new_vertices{4});
            new_edges{12} = self.create_edge(new_vertices{1}, new_vertices{5});
            
            
            %Lower left subcell
            vertices = {cll.vertices{1}, new_vertices{2}, new_vertices{1},...
                new_vertices{5}};
            edges = {new_edges{1}, new_edges{9}, new_edges{12}, new_edges{8}};
            self.create_cell(vertices, edges);
            
            %Lower right subcell
            vertices = {new_vertices{2}, cll.vertices{2}, new_vertices{3},...
                new_vertices{1}};
            edges = {new_edges{2}, new_edges{3}, new_edges{10}, new_edges{9}};
            self.create_cell(vertices, edges);
            
            %Upper right subcell
            vertices = {new_vertices{1}, new_vertices{3}, cll.vertices{3},...
                new_vertices{4}};
            edges = {new_edges{10}, new_edges{4}, new_edges{5}, new_edges{11}};
            self.create_cell(vertices, edges);
            
            %Upper left subcell
            vertices = {new_vertices{5}, new_vertices{1}, new_vertices{4},...
                cll.vertices{4}};
            edges = {new_edges{12}, new_edges{11}, new_edges{6}, new_edges{7}};
            self.create_cell(vertices, edges);
            
            %Delete old edges and cell
%             cll.delete_edges(1:4);
            delete(cll);
            
            %Update statistics
            self.nEdges = self.nEdges - 4;
            self.nCells = self.nCells - 1;
        end
                
        function map = map2fine(self, toMesh)
            %Computes map from rectangular mesh to finer regular rectangular
            %mesh given by grid vectors gridX, gridY
            %For use to map from random field discretization to PDE discret.

            %Loop through fine (FEM) cells specified by grid vectors and check
            %which coarse cell they are in
            map = zeros(toMesh.nCells, self.nCells);
            n = 1;
            for cll = self.cells
                if isvalid(cll{1})
                    m = 1;  %toMesh cell index
                    for to_cll = toMesh.cells
                        if isvalid(to_cll{1})
                            isin = cll{1}.inside(to_cll{1}.centroid);
                            if isin
                                %FEM cell m is in random fiend cell n
                                map(m, n) = 1;
                            end
                            m = m + 1;
                        end
                    end
                    n = n + 1;
                end
            end
        end
        
        function [ind, nrow, ncol] = indexIndicator(self, resolution)
            %gives boolean indices for pixels in cell. ONLY FOR SQUARE DOMAIN!
            
            if nargin < 2
                resolution = 256;
            end
            
            x = linspace(0, 1, resolution + 1);
            x = .5*(x(2:end) + x(1:(end - 1)));
            [xx, yy] = meshgrid(x);
            
            m = 1;
            for cll = self.cells
                if isvalid(cll{1})
                    ind{m} = reshape(cll{1}.inside([xx(:), yy(:)]), ...
                        resolution, resolution);
                    if nargout > 1
                        ncol(m) = max(sum(ind{m}, 2));
                        nrow(m) = max(sum(ind{m}, 1));
                    end
                    m = m + 1;
                end
            end
        end
        
                %% depreceated
        function map = map2fine_old(self, gridX, gridY)
            %Computes map from rectangular mesh to finer regular rectangular 
            %mesh given by grid vectors gridX, gridY
            %For use to map from random field discretization to PDE discret.
            
            if nargin == 2
                gridY = gridX;
            end
            
            nx = numel(gridX); ny = numel(gridY);
            cumsumX = cumsum([0, gridX]); cumsumY = cumsum([0, gridY]);
            
            
            %Loop through fine (FEM) cells specified by grid vectors and check
            %which coarse cell they are in
            map = zeros(nx*ny, self.nCells);
            n = 1;
            for cll = self.cells
                if isvalid(cll{1})
                    m = 1;  %fine (FEM) cell index
                    for y = 1:ny
                        for x = 1:nx
                            %Compute centroid of fine (FEM) cell
                            vertex_low = [cumsumX(x), cumsumY(y)];
                            vertex_high = [cumsumX(x + 1), cumsumY(y + 1)];
                            centroid = .5*(vertex_low + vertex_high);
                            isin = cll{1}.inside(centroid);
                            if isin
                                %FEM cell m is in random fiend cell n
                                map(m, n) = 1;
                            end
                            m = m + 1;
                        end
                    end
                    n = n + 1;
                end
            end
        end
    end
end

