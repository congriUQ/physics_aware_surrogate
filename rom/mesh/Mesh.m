classdef Mesh < handle
    %Mesh parent class
    properties(SetAccess = private)
        vertices            %cell array of vertex objects.
                            %Must be sorted according to global vertex number!
        edges               %cell array of edge objects.
                            %Must be sorted according to global edge number!
        cells               %cell array of cell objects.
                            %Must be sorted according to global cell number!
                            
        nVertices = 0       %Total number of vertices
    end
    
    properties(SetAccess = protected)
        nEdges = 0          %total number of edges
        nCells = 0          %total number of cells
    end
    
    methods
        function vtx = create_vertex(self, coordinates)
            %Creates a vertex and appends it to vertex list
            vtx = Vertex(coordinates);
            self.vertices{end + 1} = vtx;
            self.nVertices = self.nVertices + 1;
        end
        
        function edg = create_edge(self, vertex1, vertex2)
            %Creates an edge between vertices specified by vertex indices
            %and adds it to edge list
            edg = Edge(vertex1, vertex2);
            self.edges{end + 1} = edg;
            
            %Link vertices to edge
            vertex1.add_edge(edg);
            vertex2.add_edge(edg);
            
            self.nEdges = self.nEdges + 1;
        end
        
        function create_cell(self, vertices, edges)
            %Creates a cell and adds it to cell list
            
            self.cells{end + 1} = Cell(vertices, edges);
            
            %Link vertices to cell
            for n = 1:numel(vertices)
                vertices{n}.add_cell(self.cells{end});
            end
            
            %Link edges to cell
            for n = 1:numel(edges)
                edges{n}.add_cell(self.cells{end});
            end
            
            self.nCells = self.nCells + 1;
        end
                
        function [pv, pe, pc, fig] = plotMesh(self, fig)
            if nargin < 2
                fig = figure;
            end
            
            %plot vertices
            pv = zeros(1, numel(self.vertices));
            for n = 1:numel(self.vertices)
                hold on;
                if isvalid(self.vertices{n})
                    pv(n) = self.vertices{n}.plot(fig);
                end
            end
            
            %plot edges
            pe = zeros(1, numel(self.edges));
            for n = 1:numel(self.edges)
                hold on;
                if isvalid(self.edges{n})
                    pe(n) = self.edges{n}.plot(fig);
                end
            end
            
            %plot cells
            pc = zeros(1, numel(self.cells));
            for n = 1:numel(self.cells)
                hold on;
                if isvalid(self.cells{n})
                    c = self.cells{n}.centroid;
                    pc(n) = text(c(1), c(2), num2str(n));
                end
            end
        end
    end
end

