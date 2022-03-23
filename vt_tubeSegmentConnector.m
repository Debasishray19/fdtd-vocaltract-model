function [PV_N] = vt_tubeSegmentConnector(PV_N, prevGetRadius, currGetRadius, gridCellTypes, tubeX)
    % Store simulation grid size
    gridSize = size(PV_N);
            
    % Move through the current grid planes and store the grid cells'
    % property in currGridPlaneProp
    [currGridPlaneProp, ~] = vt_findCellTypes(PV_N, gridCellTypes, tubeX);
    
    % Move through the previous grid planes and store the grid cells'
    % property in prevGridPlaneProp
    [prevGridPlaneProp, gridCellTypeInplane] = vt_findCellTypes(PV_N, gridCellTypes, tubeX-1);
        
    % Based on the adjacent grid planes condition connect them
    if currGetRadius > prevGetRadius
        for yCount = 1:gridSize(1)
            for zCount = 1:gridSize(3)
                if (currGridPlaneProp(yCount, zCount) == gridCellTypeInplane.inVTContour && prevGridPlaneProp(yCount, zCount) == gridCellTypeInplane.outVTContour)
                    PV_N(yCount, tubeX-1, zCount, 5) = gridCellTypes.cell_wall;              
                end
            end 
        end
    else
       for yCount = 1:gridSize(1)
            for zCount = 1:gridSize(3)
                if (prevGridPlaneProp(yCount, zCount) == gridCellTypeInplane.inVTContour && currGridPlaneProp(yCount, zCount) == gridCellTypeInplane.outVTContour)
                    PV_N(yCount, tubeX, zCount, 5) = gridCellTypes.cell_wall;               
                end
            end 
       end         
    end
end