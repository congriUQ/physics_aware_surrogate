classdef Edge < handle
    %Mesh edge class
    
    properties(SetAccess = private)
        vertices                %adjacent vertices
        cells                   %adjacent cells
        length                  %length of edge element
    end
    
    methods
        function self = Edge(vertex1, vertex2)
            %Constructor
            self.vertices = {vertex1, vertex2};
            self.length = norm(vertex1.coordinates - vertex2.coordinates);
        end
        
        function self = add_cell(self, cell)
            self.cells{end + 1} = cell;
        end
        
        function [l] = plot(self, fig)
            if nargin > 1
                figure(fig)
            end
            
            hold on;
            c1 = self.vertices{1}.coordinates;
            c2 = self.vertices{2}.coordinates;
            l= line([c1(1) c2(1)], [c1(2) c2(2)]);
        end
    end
end

