function [shortestPath] = shortestPath(lambdak, dir, resolution, distMeasure)
%Gives minimal distance of fluid phase path from left to right (dir == 'x')
%or top to bottom (dir == 'y')
%   lambdak:        binary pixel input. true is exclusion, false is fluid
%   distMeasure:    'cityblock' or 'quasi-euclidean', 'chessboard' gives
%                   the edge length always!
if nargin < 3
    resolution = 256;
end

%fluid is true, exclusion is false
lambdak = ~lambdak;

shortestPath = Inf;
if(dir == 'x')
    %Check first if there are true elements on the relevant boundaries
    %(left/right for x)
    left = lambdak(:, 1);
    right = lambdak(:, end);
    if(any(left) && any(right))
        %Loop through left boundary and use right boundary as mask for
        %matlab bwgeodesic function
%         leftIndex = find(left);
        %only consider middle pixel -- otherwise this is always == edge length
        %for 2x2 model
        leftIndex = round(size(left, 1)/2);
        incr = 0;
        while(~isfinite(shortestPath) && leftIndex > 1 &&...
                leftIndex < size(left, 1))
            leftIndex = leftIndex + incr;
            incr = -(incr + 1);
            geo = bwdistgeodesic(lambdak, 1, leftIndex, distMeasure);
            if(any(isfinite(geo(:, end))))
                shortestPathTemp = min(geo(:, end));
                if(shortestPathTemp < shortestPath)
                    shortestPath = shortestPathTemp;
                    if(~isfinite(shortestPath))
                        leftIndex
                        dir
                        geo
                        error('Zero path length of connected path')
                    end
                end
            end
        end
        if ~isfinite(shortestPath)
            shortestPath = pi*size(left, 1);
        end
    end
elseif(dir == 'y')
    %Check first if there are true elements on the relevant boundaries 
    %(top/bottom for y)
    top = lambdak(1, :);
    bottom = lambdak(end, :);
    if(any(top) && any(bottom))
        %Loop through upper boundary and use lower boundary as mask for matlab
        %bwgeodesic function
%         topIndex = find(top);
        %only consider middle pixel -- otherwise this is always == edge length
        %for 2x2 model
        topIndex = round(size(top, 2)/2);
        incr = 0;
        while(~isfinite(shortestPath) && topIndex > 1 &&...
                topIndex < size(top, 2))
            topIndex = topIndex + incr;
            incr = -(incr + 1);
            geo = bwdistgeodesic(lambdak, topIndex, 1, distMeasure);
            if(any(isfinite(geo(end, :))))
                shortestPathTemp = min(geo(end, :));
                if(shortestPathTemp < shortestPath)
                    shortestPath = shortestPathTemp;
                    if(~isfinite(shortestPath))
                        topIndex
                        dir
                        geo
                        error('Zero path length of connected path')
                    end
                end
            end
        end
        if ~isfinite(shortestPath)
            shortestPath = pi*size(top, 2);
        end
    end
else
    error('which direction?')
end
%the +1 is heuristic -- why is the shortest path sometimes one pixel
%smaller than cell edge length?
shortestPath = (shortestPath + 1)/resolution;
end

