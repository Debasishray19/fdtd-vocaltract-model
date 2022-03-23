function [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo] = vt_ellipticalTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds)
 
    % Define units
    meter = 1;
    centimeter = 1e-2*meter;
    milimeter = 1e-3*meter;
    
    % STEP XX: Define tube geometry
    [segment_length, tubeSectionArea_incm2] = vt_areaFunction(simulationType, vowel);
 
    % Tube section area in m^2
    tubeSectionArea_inm2 = tubeSectionArea_incm2.*(centimeter*centimeter);
    
    % Virtual mic/listener position near the mouth-end [inside the tube]
    micXdistanceFromMouthEnd = 3*milimeter;
    micxdistanceInCellsFromMouthEnd = round(micXdistanceFromMouthEnd/ds);

    % number of tube segments
    numSections = length(tubeSectionArea_inm2);

    % Calculate tube cross-section semi-major and semi-minor axis length
    % The length of tube major and minor axis should be computed from the
    % original tube cross-sections/geometry [Arnela et al. 2016 - Influence 
    % of vocal tract geometry simplifications]. However, in the case of
    % missing geometry images, the length of the semi-major and semi-minor 
    % axis can be determined by defining a ratio by the user.
    
    isTubeGeometryDataExist = 0;
    
    if isTubeGeometryDataExist
        % TODO: We need vocal tract tube geometry/MRI images for this
    else
        % [Note: Keep semimajorAxisRatio >  semiminorAxisRatio]
        semimajorAxisRatio = 3;
        semiminorAxisRatio = 1;
        fprintf("Ratio of the semi-major:semi-minor = %d:%d\n", semimajorAxisRatio, semiminorAxisRatio);
        
        % Determine length of the semi-major and semi-minor axis length
        % if semimajorAxis:semiminorAxis = a:b then their length can be
        % semimajorAxisLen = ax,  semiminorAxisLen = bx
        % ellipseArea = pi*ax*bx;
        % x = sqrt(ellipseArea/(semimajorAxisLen*semiminorAxisLen))
        
        lenX = sqrt(tubeSectionArea_inm2./(semimajorAxisRatio*semiminorAxisRatio*pi));
        majorAxisLen = lenX.*(semimajorAxisRatio*2);
        minorAxisLen = lenX.*(semiminorAxisRatio*2);
    end
    
    tubeMajorAxisLenInCells = round(majorAxisLen./ds);
    tubeMinorAxisLenInCells = round(minorAxisLen./ds);
        
    % Cross-sectional area at the tube start end
    % I didn't use "tubeSectionArea_inm2" directly because I wanted to use
    % the derived area
    tubeStartArea = pi*(majorAxisLen(1)/2)*(minorAxisLen(1)/2);
    
    % Change the majorAxisLenInCells and minorAxisLenInCells to 1 
    % if it contains 0
    tubeMajorAxisLenInCells(tubeMajorAxisLenInCells==0)=1;
    tubeMinorAxisLenInCells(tubeMinorAxisLenInCells==0)=1;
     
    % STEP XX: Choose the best possible odd number from the MinorAxisLenInCell 
    % and tubeMajorAxisLenInCells array. In this model, the major axis goes 
    % along the xz-plane.
    for counter = 1:numSections
        
        % Verify if the tubeMajorAxisLenInCells is odd or not
        if mod(tubeMajorAxisLenInCells(counter), 2) == 0
            
            % Find the difference between the rounded and the actual
            % diameter value = Estimated diameter-Actual diameter
            diff = tubeMajorAxisLenInCells(counter) - ...
                (majorAxisLen(counter)/ds);
            
            if diff>0
                tubeMajorAxisLenInCells(counter) = ...
                    tubeMajorAxisLenInCells(counter)-1;
            else
                tubeMajorAxisLenInCells(counter) = ...
                    tubeMajorAxisLenInCells(counter)+1;
            end    
        end
        
        % Verify if the tubeMinorAxisLenInCells is odd or not
        if mod(tubeMinorAxisLenInCells(counter), 2) == 0
            
            % Find the difference between the rounded and the actual
            % diameter value = Estimated diameter-Actual diameter
            diff = tubeMinorAxisLenInCells(counter) - ...
                (minorAxisLen(counter)/ds);
            
            if diff>0
                tubeMinorAxisLenInCells(counter) = ...
                    tubeMinorAxisLenInCells(counter)-1;
            else
                tubeMinorAxisLenInCells(counter) = ...
                    tubeMinorAxisLenInCells(counter)+1;
            end    
        end        
    end
    
    % Store the ellipse semi-major axis and semi-minor axis diameter
    ellipseAxisLenInfo = zeros(2, numSections);
    ellipseAxisLenInfo(1,:) = tubeMajorAxisLenInCells;
    ellipseAxisLenInfo(2,:) = tubeMinorAxisLenInCells;
    
    % STEP XX: Find the total tube length and calculate the percentage error 
    % in the approximated tube length

    % Number of cells for total tube length
    actualTubeLength = numSections*segment_length;
    totalTubeLengthInCells = round(actualTubeLength/ds);
    approxTubeLength = totalTubeLengthInCells*ds;
    
    % Percentage error in total tube length
    totalTubeLengthError = (approxTubeLength-actualTubeLength)/...
                           (actualTubeLength);
    fprintf('Approximated tube length percentage error = %.4f \n',totalTubeLengthError*100);
        
    if baffleSwitchFlag == 1
        % TO DO: Deb Implement the baffle condition
        % Deb: For eccentric tube condition how should we include baffle?
    else
        % The major axis of the elliptical cross-section goes along the
        % z-axis and the minor axis of the elliptical cross-section goes
        % along the y-axis.
        domainX = totalTubeLengthInCells + 1 + 1; % +1 is excitation and and dirichlet layer
        domainY = max(tubeMinorAxisLenInCells) + 2; % +2 is tube walls
        domainZ = max(tubeMajorAxisLenInCells) + 2; % +2 is tube walls
    end
    
    % Derive frame size
    frameX = domainX + 2; %2 = dead cells
    frameX = frameX + 2*pmlSwitch*pmlLayer;
    
    frameY = domainY + 2; %2 = dead cells
    frameY = frameY + 2*pmlSwitch*pmlLayer;
    
    frameZ = domainZ + 2; %2 = dead cells
    frameZ = frameZ + 2*pmlSwitch*pmlLayer;
    
    % STEP XX: Add all the frame params to the simGridParam
    simGridParam.frameX = frameX;
    simGridParam.frameY = frameY;
    simGridParam.frameZ = frameZ;
    
    % Define PV_N array to store pressure and velocity components
    % PV_N(:,:,1) = To store cell pressure
    % PV_N(:,:,2) = To store Vx
    % PV_N(:,:,3) = To store Vy
    % PV_N(:,:,4) = To store Vz
    % PV_N(:,:,5) = To store grid cell types
    PV_N = zeros(frameY, frameX, frameZ, 5);
    
    % Define cell types and store it in PV_N(,,5)
    % Declare all the cells as air by default
    PV_N(1:frameY, 1:frameX, 1:frameZ, 5) = gridCellTypes.cell_air;
    
    % STEP XX: Create the regular tube geometry
    
    % The axis of the elliptical tube is along the x-axis, direction along
    % which acoustic wave propagates. We will start building the tube model
    % from the starting position towards the end position. 
    
    % STEP XX.XX Find the center of the tube along the yz-plane
    centerY = ceil(frameY/2);
    centerZ = ceil(frameZ/2);
    
    % STEP XX.XX Find the starting and ending point of the tube.
    % 1 is: for starting point, dead layer and excitation
    % 2 is: for cell_head in case of the circular baffle and to have an 
    % empty cell between the cell_head and pml layers.
    if rad ~=3
        tubeStart.startX = 1 + 1 + pmlLayer*pmlSwitch + 1;
        tubeEnd.endX     = tubeStart.startX + totalTubeLengthInCells - 1 ;
    else
        %TO DO: Deb check this later tubeStartX = 1 + 1 + pmlLayer*pmlSwitch + 1+2;
    end
    tubeStart.startY = centerY;
    tubeEnd.endY     = centerY;
    
    tubeStart.startZ = centerZ;
    tubeEnd.endZ     = centerZ;
    
    % STEP XX: Store the cummulative length of each tube section
    tubeCummSectionLength = zeros(1, numSections);
    for sectionCount = 1:numSections
        tubeCummSectionLength(sectionCount) = segment_length*sectionCount;
    end
    
    [PV_N, currTubeSectionDiameterCells_SegmentCounter]...
        = vt_ellipticalContourGeneration(PV_N, numSections, totalTubeLengthInCells, tubeCummSectionLength, ellipseAxisLenInfo, boundaryInterpolation, tubeStart, gridCellTypes, ds);
       
    % STEPXX: Define excitation cells
    % Check the grid cell types for the yz-plane that exists besides
    % (left-side) the tube starting point
    xExcitationPos = tubeStart.startX-1;
    [gridPlaneProp, gridCellTypeInplane] = vt_findCellTypes(PV_N, gridCellTypes, tubeStart.startX);
    
    % Find the grid size and traverse through yz-plane to assign
    % cell_excitation
    gridSize = size(PV_N);
    
    for yCount = 1:gridSize(1)
        for zCount = 1:gridSize(3)
            if gridPlaneProp(yCount, zCount) == gridCellTypeInplane.inVTContour
                PV_N(yCount, xExcitationPos, zCount, 5) = gridCellTypes.cell_excitation;
            elseif gridPlaneProp(yCount, zCount) == gridCellTypeInplane.onVTContour
                PV_N(yCount, xExcitationPos, zCount, 5) = gridCellTypes.cell_wall;
            end
        end
    end
    
    % STEPXX: Define no_Pressure cells [Dirichlet Boundary Condition]   
    if rad==2
        tubeEnd.endX = tubeStart.startX + totalTubeLengthInCells - 1 ;
        [gridPlaneProp, gridCellTypeInplane] = vt_findCellTypes(PV_N, gridCellTypes, tubeEnd.endX);
        
        % Traverse through yz-plane to assign cell_noPressure
        for yCount = 1:gridSize(1)
            for zCount = 1:gridSize(3)
                if gridPlaneProp(yCount, zCount) == gridCellTypeInplane.inVTContour
                    PV_N(yCount, tubeEnd.endX+1, zCount, 5) = gridCellTypes.cell_noPressure;
                elseif gridPlaneProp(yCount, zCount) == gridCellTypeInplane.onVTContour
                    PV_N(yCount, tubeEnd.endX+1, zCount, 5) = gridCellTypes.cell_noPressure;
                end
            end
        end
    end
    
    % Place two microphines:
    % 1. Near mouth-end [3mm inside of mouth] = listenerInfo
    % 2. Near glottal-end [Next to excitation cells] = sourceInfo
    
    % STEPXX: Determine listener position   
    listenerInfo.listenerX = tubeEnd.endX-micxdistanceInCellsFromMouthEnd;
    listenerInfo.listenerY = tubeEnd.endY;
    listenerInfo.listenerZ = tubeEnd.endZ; 
    
    % STEPXX: Determine the source position near the glottal-end  
    sourceInfo.sourceX = tubeStart.startX+1;
    sourceInfo.sourceY = tubeStart.startY;
    sourceInfo.sourceZ = tubeStart.startZ;
end