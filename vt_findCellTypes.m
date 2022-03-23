function [gridPlaneProp, gridCellTypeInplane] = vt_findCellTypes(PV_N, gridCellTypes, tubeX)

    % For a given yz-plane, this function determines whether a particular
    % grid cell is inside the tube contour, outside of the tube contour or 
    % on the tube contour [tubeContour].
    
    % Define grid cells as inside the VT contour, outside the VT contour or
    % on the VT contour.
    gridCellTypeInplane.inVTContour = 1;
    gridCellTypeInplane.outVTContour = 2;
    gridCellTypeInplane.onVTContour = 3;
    
    gridSize = size(PV_N);
    gridPlaneProp = zeros(gridSize(1), gridSize(3));
    
    for zCount = 1:gridSize(3)
        findPlaneWalls = find(PV_N(:,tubeX, zCount, 5) == gridCellTypes.cell_wall);
        gridPlaneProp(findPlaneWalls, zCount) = gridCellTypeInplane.onVTContour;
        if isempty(findPlaneWalls)
            gridPlaneProp(:, zCount) = gridCellTypeInplane.outVTContour;
        else
            minWall = min(findPlaneWalls);
            maxWall = max(findPlaneWalls);
            gridPlaneProp(1:minWall-1, zCount) = gridCellTypeInplane.outVTContour;
            gridPlaneProp(maxWall+1:gridSize(1), zCount) = gridCellTypeInplane.outVTContour;
            
            inVTContourCells = gridPlaneProp(:, zCount) == 0;
            gridPlaneProp(inVTContourCells, zCount) = gridCellTypeInplane.inVTContour;
        end
    end
end