function  [PV_N, currTubeSectionDiameterCells_SegmentCounter]  = vt_eccentricEllipticalContourGeneration(PV_N, numSections, totalTubeLengthInCells, tubeCummSectionLength, ellipseAxisLenInfo, boundaryInterpolation, tubeStart, pmlSwitch, pmlLayer, gridCellTypes, ds)

    % Set a counter to traverse through tubeCummSectionLength
    tubeSegmentCounter = 1;
    cellHalfLen = ds/2;
    
    % Define currTubeSectionDiameterCells
    % For ellptical cross-section, we will use this array to store the minor 
    % axis tube diamaters. 
    currTubeSectionDiameterCells_SegmentCounter = zeros(2, totalTubeLengthInCells);
    tubeSemiMinorAxisRadiusInCells = zeros(1, totalTubeLengthInCells);
    tubeSemiMajorAxisRadiusInCells = zeros(1, totalTubeLengthInCells);
    
    % Set starting coordinte of tube mid-sagittal axis
    startX = tubeStart.startX;
    startZ = tubeStart.startZ;
    
    % STEP1: Each tube segment consists of number of yz-planes. We first store the
    % tube semi-minor axis radius/length for each of those planes.
    
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
            currSemiMajorAxisRadius = (ellipseAxisLenInfo(1, tubeSegmentCounter)-1)/2;
            currSemiMinorAxisRadius = (ellipseAxisLenInfo(2, tubeSegmentCounter)-1)/2;
            
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
                currSemiMajorAxisRadius = (ellipseAxisLenInfo(1, tubeSegmentCounter)-1)/2;
                currSemiMinorAxisRadius = (ellipseAxisLenInfo(2, tubeSegmentCounter)-1)/2;
                
                % store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
            else
                % Get the radius for that cross-section [-1: each segment consist of odd number of cells]
                currSemiMajorAxisRadius = (ellipseAxisLenInfo(1, tubeSegmentCounter)-1)/2;
                currSemiMinorAxisRadius = (ellipseAxisLenInfo(2, tubeSegmentCounter)-1)/2;
                
                % store the current cross-section tube segment counter
                currTubeSectionDiameterCells_SegmentCounter(2,tubeLenCellsCount) = tubeSegmentCounter;
        
                % And then increase the tubeSegmentCounter 
                tubeSegmentCounter = tubeSegmentCounter+1;
            end
        end       
        
        tubeSemiMinorAxisRadiusInCells(tubeLenCellsCount) = currSemiMinorAxisRadius;
        tubeSemiMajorAxisRadiusInCells(tubeLenCellsCount) = currSemiMajorAxisRadius;
    end
    
    if boundaryInterpolation == 1
        % TODO: Deb - Think of what's the correct method to interpolate
        % between two tube segmnts.
    end
    
    % Store the current cross-section tube segment semi-minor axis diameter
    currTubeSectionDiameterCells_SegmentCounter(1,:) = tubeSemiMinorAxisRadiusInCells.*2+1;
    
    % STEP2: Now as we have tubeRadius for each yz-plane, draw the vocal tract
    % contour starting from the tubeStart.startX
    
    for tubeLenCellsCount = 1:totalTubeLengthInCells
        currSemiMajorAxisRadius = tubeSemiMajorAxisRadiusInCells(tubeLenCellsCount);
        currSemiMinorAxisRadius = tubeSemiMinorAxisRadiusInCells(tubeLenCellsCount);
        
        % The Y and Z coordinate of the tube central axis is not gonna
        % change. Only the X axis will change as we move away from the 
        % glottal-end. We need to draw an ellipse in the yz-plane
        
        % [-1 to start constructing tube exactly from the tubeStart.startX
        % position]
       
        % Deb: You mean tubeUpperY and tubeLowerY in yz-plane
        tubeX = startX + (tubeLenCellsCount-1); %[For each iteration tubeX remains fix]
        
        tubeUpperY = 1+(pmlLayer*pmlSwitch)+1;  % deadcell+pmllayers +tube Wall
        tubeLowerY = tubeUpperY+(currSemiMinorAxisRadius*2+1)+1; % tubeUpperY+tubeDiameter+nextcell
        tubeMidY   = tubeUpperY+currSemiMinorAxisRadius+1;
        
        % Deb: tubeRightZ and tubeLeftZ mean the right and left side of
        % the center in yz-plane.
        tubeRightZ = startZ;
        tubeLeftZ = startZ;
        
        PV_N(tubeUpperY, tubeX, tubeRightZ, 5) = gridCellTypes.cell_wall;
        PV_N(tubeLowerY, tubeX, tubeLeftZ, 5) = gridCellTypes.cell_wall;
        
        % Store the yCoordinate to close if there are any gaps
        prevtubeUpperY = tubeUpperY;
        prevtubeLowerY = tubeLowerY;
        
        for semiMajorAxisRadiusCounter = 1:currSemiMajorAxisRadius
            % Ellipse equation: z^2/a^2 + y^2/b^2 = 1 [Note: I have mentioned 
            % z^2 becasuse the semi-major axis is along the z-axis]
            % a = (currSemiMajorAxisRadius+0.5)*ds [Length of semiMaorAxisRadius]
            % b = (currSemiMinorAxisRadius+0.5)*ds [Length of semiMinorAxisRadius]
            % z = semiMajorAxisRadiusCounter*ds
            % Determine y (ellipse_height) = ? using Ellipse equation [i.e., height of the ellipse along y-axis]
            
            a = (currSemiMajorAxisRadius+0.5)*ds;
            b = (currSemiMinorAxisRadius+0.5)*ds;
            z = semiMajorAxisRadiusCounter*ds;
            
            currEllipseHeightFromCentralAxis = b*sqrt(1- (z^2/a^2));
            currHeightFromTopOfCentralAxisRowCell = currEllipseHeightFromCentralAxis-(ds/2);
            
            currHeightInCells = round(currHeightFromTopOfCentralAxisRowCell/ds);
            
            % Define tube cross-section coordinates in the yz plane
            tubeRightZ = startZ+semiMajorAxisRadiusCounter;
            tubeLeftZ  = startZ-semiMajorAxisRadiusCounter;
            
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
        tubeRightZ = startZ+semiMajorAxisRadiusCounter+1;
        tubeLeftZ  = startZ-semiMajorAxisRadiusCounter-1;
        
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
        if tubeX ~= startX && prevSemiMajorAxisRadius   ~= currSemiMajorAxisRadius          
            PV_N = vt_tubeSegmentConnector(PV_N, prevSemiMajorAxisRadius, currSemiMajorAxisRadius, gridCellTypes, tubeX);
        end
        
        % Store the gridRadius for the current yz-plane
        prevSemiMajorAxisRadius = currSemiMajorAxisRadius;
    end   
end