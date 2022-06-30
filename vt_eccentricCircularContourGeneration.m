function  [PV_N, currTubeSectionDiameterCells_SegmentCounter]  = vt_eccentricCircularContourGeneration(PV_N, numSections, totalTubeLengthInCells, tubeSectionDiameterCells, tubeCummSectionLength, boundaryInterpolation, tubeStart, pmlSwitch, pmlLayer, gridCellTypes, ds)

    % Set a counter to traverse through tubeCummSectionLength
    tubeSegmentCounter = 1;
    cellHalfLen = ds/2;
    
    % Define currTubeSectionDiameterCells
    currTubeSectionDiameterCells_SegmentCounter = zeros(2, totalTubeLengthInCells);
    tubeRadiusArray = zeros(1, totalTubeLengthInCells);
    
    % Set starting coordinte of tube mid-sagittal axis
    % For eccentric geometry we don't need startY
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
            currGetRadius = (tubeSectionDiameterCells(tubeSegmentCounter)-1)/2;
            
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
                currGetRadius = (tubeSectionDiameterCells(tubeSegmentCounter)-1)/2;
                
                % store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
            else
                % Get the radius for that cross-section [-1: each segment consist of odd number of cells]
                currGetRadius = (tubeSectionDiameterCells(tubeSegmentCounter)-1)/2;
                
                % store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
        
                % And then increase the tubeSegmentCounter 
                tubeSegmentCounter = tubeSegmentCounter+1;
            end
        end       
        
        tubeRadiusArray(tubeLenCellsCount) = currGetRadius;
    end
    
    if boundaryInterpolation == 1
        [tubeRadiusArray] = vt_boundaryInterpolation(tubeRadiusArray, ds);
    end
    
    % Store the current cross-section tube segment diameter
    currTubeSectionDiameterCells_SegmentCounter(1,:) = tubeRadiusArray.*2+1;
    
    % STEP2: Now as we have tubeRadius for each yz-plane, draw the vocal tract
    % contour starting from the tubeStart.startX
    
    for tubeLenCellsCount = 1:totalTubeLengthInCells
        currGetRadius = tubeRadiusArray(tubeLenCellsCount);
        
        % For eccentric configuration, the y coordinate of the tube central
        % axis will change but the z coordinate will remain the same for
        % all the tube segments
        % And the X axis will change as we move away from the 
        % glottal-end. We need to draw a circle in the yz-plane
        
        % [-1 to start constructing tube exactly from the tubeStart.startX
        % position]
       
        % Deb: You mean tubeUpperY and tubeLowerY in yz-plane
        tubeX = startX + (tubeLenCellsCount-1); %[For each iteration tubeX remains fix]
        
        tubeUpperY = 1+(pmlLayer*pmlSwitch)+1;  % deadcell+pmllayers +tube Wall
        tubeLowerY = tubeUpperY+(currGetRadius*2+1)+1; % tubeUpperY+tubeDiameter+nextcell
        tubeMidY   = tubeUpperY+currGetRadius+1;
        
        % Deb: tubeRightZ and tubeLeftZ mean the right and left side of
        % the center in yz-plane.
        tubeRightZ = startZ;
        tubeLeftZ = startZ;
  
        PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = gridCellTypes.cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 5) = gridCellTypes.cell_wall;
        
        % Store the yCoordinate to close if there are any gaps
        prevtubeUpperY = tubeUpperY;
        prevtubeLowerY = tubeLowerY;
        
        for radiusCounter = 1:currGetRadius
            % rActual = (currGetRadius+0.5)ds
            % zPrime = (radiusCounter)ds

            rActual = (currGetRadius+0.5)*ds;
            zPrime = (radiusCounter)*ds;
            currHeightFromCentralAxis = sqrt(rActual^2 - zPrime^2);
            currHeightFromTopOfCentralAxisRowCell = currHeightFromCentralAxis-(ds/2);
            
            currHeightInCells = round(currHeightFromTopOfCentralAxisRowCell/ds);
            
            % Define tube cross-section coordinates in the yz plane
            tubeRightZ = startZ+radiusCounter;
            tubeLeftZ  = startZ-radiusCounter;
            
            % Define tubeUpperY and tubeLowerY
            tubeUpperY = tubeMidY-currHeightInCells-1;
            tubeLowerY = tubeMidY+currHeightInCells+1;
            
            PV_N(tubeUpperY,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
            PV_N(tubeLowerY,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
            
            PV_N(tubeUpperY,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
            PV_N(tubeLowerY,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
            
            % Connect grid cells in yz-plane if there are any gaps
            if tubeUpperY-prevtubeUpperY > 1
                PV_N(prevtubeUpperY+1:tubeUpperY,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
                PV_N(prevtubeUpperY+1:tubeUpperY,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
                PV_N(tubeLowerY:prevtubeLowerY-1,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
                PV_N(tubeLowerY:prevtubeLowerY-1,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
            end
            
            prevtubeUpperY = tubeUpperY;
            prevtubeLowerY = tubeLowerY;      
        end
        
        % In the end close the circle by using walls at the complete right
        % and left side of the yz-plane [as per the tube radius]
        
        tubeUpperY = tubeMidY;
        tubeLowerY = tubeMidY;
        tubeRightZ = startZ+currGetRadius+1;
        tubeLeftZ  = startZ-currGetRadius-1;
        
        PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = gridCellTypes.cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 5)  = gridCellTypes.cell_wall;
        
        if tubeUpperY-prevtubeUpperY > 1
            PV_N(prevtubeUpperY+1:tubeUpperY,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
            PV_N(prevtubeUpperY+1:tubeUpperY,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
            PV_N(tubeLowerY:prevtubeLowerY-1,tubeX,tubeRightZ,5) = gridCellTypes.cell_wall;
            PV_N(tubeLowerY:prevtubeLowerY-1,tubeX,tubeLeftZ,5) = gridCellTypes.cell_wall;
        end
        
        % Connect grid cells between two adjacent yz-planes if there are 
        % any gaps
        if tubeX ~= startX && prevGetRadius ~= currGetRadius          
            PV_N = vt_tubeSegmentConnector(PV_N, prevGetRadius, currGetRadius, gridCellTypes, tubeX);
        end
        
        % Store the gridRadius for the current yz-plane
        prevGetRadius = currGetRadius;
    end   
end