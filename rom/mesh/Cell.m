classdef Cell < handle
    %Mesh cell class
    
    properties(SetAccess = private)
        vertices                %adjacent vertices. Have to be sorted according
                                %to local vertex number!(counterclockwise
                                %starting from lower left corner)
        edges                   %adjacent edges. Have to be sorted according to
                                %local edge number! (counterclockwise starting 
                                %from lower left corner)
        surface                 %surface of cell element
        centroid                %centroid of cell
        type = 'rectangle'
    end
    
    methods
        function self = Cell(vertices, edges, type)
            %Constructor
            %Vertices and edges must be sorted according to local
            %vertex/edge number!
            self.vertices = vertices;
            self.edges = edges;
            if nargin > 2
                self.type = type;
            end
            
            %Compute centroid
            self.centroid = 0;
            for n = 1:numel(self.vertices)
                self.centroid = self.centroid + self.vertices{n}.coordinates;
            end
            self.centroid = self.centroid/numel(self.vertices);
            
            self.compute_surface;
        end
        
        function delete_edges(self, indices)
            %Deletes edges according to indices
            for i = indices
                delete(self.edges{i});
                self.edges{i} = [];
            end
        end
        
        function compute_surface(self)
            if strcmp(self.type, 'rectangle')
                self.surface = self.edges{1}.length*self.edges{2}.length;
            else
                error('unknown cell type')
            end
        end
        
        function [isin] = inside(self, x)
            %Check if point is inside of cell
            %   x:  N*dim array, --> each line is a point
            
            isin = (x(:, 1) > self.vertices{1}.coordinates(1) - eps & ...
                x(:, 1) <= self.vertices{3}.coordinates(1) + eps & ...
                x(:, 2) > self.vertices{1}.coordinates(2) - eps & ...
                x(:, 2) <= self.vertices{3}.coordinates(2) + eps);
            
            
            
        end
    end
end

