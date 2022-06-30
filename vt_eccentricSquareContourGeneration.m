function  [PV_N, currTubeSectionDiameterCells_SegmentCounter]  = vt_eccentricSquareContourGeneration(PV_N, numSections, totalTubeLengthInCells, tubeSquareLenInCells, tubeCummSectionLength, boundaryInterpolation, tubeStart, pmlSwitch, pmlLayer, gridCellTypes, ds)

    % Set a counter to traverse through tubeCummSectionLength
    tubeSegmentCounter = 1;
    cellHalfLen = ds/2;
    
    % Define currTubeSectionDiameterCells
    % [Note]: The variable "tubeRadiusArray" might be misleading since
    % there is no concept of radius in a square. Here the radius array
    % means height of the tube from the center.
    
    currTubeSectionDiameterCells_SegmentCounter = zeros(2, totalTubeLengthInCells);
    tubeRadiusArray = zeros(1, totalTubeLengthInCells);
    
    % Set starting coordinte of tube mid-sagittal axis
    startX = tubeStart.startX;
    startZ = tubeStart.startZ;
    
    % STEP1: Each tube segment consists of number of yz-planes. We first store the
    % tube radius for each of those planes.
    
    for tubeLenCellsCount = 1:totalTubeLengthInCells
        
        % Check the current tube length
        currTubeLength = tubeLenCellsCount*ds;
        
        % Verify if the currTubeLength is more than the tubeCummSectionLength
        % for the current sectionCounter
        % if small or equal then set the tube wall as expected-normal Case
        if currTubeLength <= tubeCummSectionLength(tubeSegmentCounter)
            
            % Get the tube Radius
            % We are subtracting 1 as we'll assume that there is a middle
            % row of cells which will act like a mirror/centerline.
            currGetRadius = (tubeSquareLenInCells(tubeSegmentCounter)-1)/2;
            
            % store the current cross-section tube segment counter
            currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
        else     
            % If the currTubeLength is greater than the actual tubeCummSectionLength
            % Find the difference between currTubeLength and tubeCummSectionLength
            diffLength = currTubeLength - tubeCummSectionLength(tubeSegmentCounter);
            
            if diffLength>cellHalfLen && tubeSegmentCounter~=numSections
                % Increase the tubeSegmentCounter
                tubeSegmentCounter = tubeSegmentCounter+1;
                
                % Get the radius for that cross-section [-1: each segment consist of odd number of cells]
                currGetRadius = (tubeSquareLenInCells(tubeSegmentCounter)-1)/2;
                
                % store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
            else
                % Get the radius for that cross-section [-1: each segment consist of odd number of cells]
                currGetRadius = (tubeSquareLenInCells(tubeSegmentCounter)-1)/2;
                
                % store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
        
                % And then increase the tubeSegmentCounter 
                tubeSegmentCounter = tubeSegmentCounter+1;
            end
        end       
        
        tubeRadiusArray(tubeLenCellsCount) = currGetRadius;
    end
    
    if boundaryInterpolation == 1
        % TODO: Deb - Think of what's the correct method to interpolate
        % between two tube segmnts.
    end
    
    % Store the current cross-section tube segment diameter
    currTubeSectionDiameterCells_SegmentCounter(1,:) = tubeRadiusArray.*2+1;
    
    % STEP2: Now as we have tubeRadius for each yz-plane, draw the vocal tract
    % contour starting from the tubeStart.startX
    
    for tubeLenCellsCount = 1:totalTubeLengthInCells
        currGetRadius = tubeRadiusArray(tubeLenCellsCount);
        
        % The Y and Z coordinate of the tube central axis is not gonna
        % change. Only the X axis will change as we move away from the 
        % glottal-end. We need to draw a square in the yz-plane
        
        % [-1 to start constructing tube exactly from the tubeStart.startX
        % position]
       
        tubeX = startX + (tubeLenCellsCount-1); %[For each iteration tubeX remains fix]
        
        % Deb: You mean tubeUpperY and tubeLowerY in yz-plane. As the tube
        % cross-section has the square shape, the y-coordinate remains
        % the same.
        tubeUpperY = 1+(pmlLayer*pmlSwitch)+1;  % deadcell+pmllayers +tube Wall
        tubeLowerY = tubeUpperY+(currGetRadius*2+1)+1; % tubeUpperY+tubeDiameter+nextcell
        
        % find the left-most and right most z-coordinate.
        tubeLeftZ = startZ-currGetRadius-1;
        tubeRightZ = startZ+currGetRadius+1;
        
        % Set the tube wall
        PV_N(tubeUpperY,tubeX,tubeLeftZ:tubeRightZ,5) = gridCellTypes.cell_wall;
        PV_N(tubeLowerY,tubeX,tubeLeftZ:tubeRightZ,5) = gridCellTypes.cell_wall;
        PV_N(tubeUpperY:tubeLowerY,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
        PV_N(tubeUpperY:tubeLowerY,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
        
        % Connect grid cells between two adjacent yz-planes if there are 
        % any gaps
        if tubeX ~= startX && prevGetRadius ~= currGetRadius         
            PV_N = vt_tubeSegmentConnector(PV_N, prevGetRadius, currGetRadius, gridCellTypes, tubeX);
        end
        
        % Store the gridRadius for the current yz-plane
        prevGetRadius = currGetRadius;
    end   
end