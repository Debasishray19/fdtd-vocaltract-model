function [tubeNewRadiusArray] = vt_boundaryInterpolation(tubeRadiusArray, ds)
    
    % Illustration:
    %                     C         C
    %                     *         *
    %             E   *   *         *   *   E
    % 		  	  *	  *   *   or    *   *   *
    % 		  *	  *	  *   *         *   *   *   *
    % 	  *	  *	  *	  *   *         *   *   *   *   *
    % *********************         ********************* 
    % A           D       B         B       D           A
    %         (Case 1)         or          (case 2)
    
    % In the Triangle ABC, let's assume, we need to find DE = ?
    % AD/AB = DE/BC => DE = (AD*BC)/AB; [Similar Triangle Equality]
    
    % Create a new tubeRadiusArray to store the radius of the interpolated
    % boundary.
    tubeNewRadiusArray = zeros(1, length(tubeRadiusArray));
    
    % Set a counter to point to the start of the triangle
    triangleStartCounter = 1;
    
    while triangleStartCounter < length(tubeRadiusArray)
        
        % Set a counter and traverse from start of the tube to the end of the tube
        tubeCounter = 1;
        
        % To use the "Similar Triangle Equality" first find the start and end of the triangle
        while tubeRadiusArray(triangleStartCounter) == tubeRadiusArray(triangleStartCounter+tubeCounter)
             if triangleStartCounter+tubeCounter ~= length(tubeRadiusArray)
                 tubeCounter = tubeCounter+1;
             else
                 break;
             end
        end
        
        % Set triangle length and height
        triangleLengthInCells = tubeCounter;
        triangleHeightInCells = tubeRadiusArray(triangleStartCounter+tubeCounter) - tubeRadiusArray(triangleStartCounter);
        
        % Set the tubeNewRadiusArray for the trainagleStart and trianglEnd
        tubeNewRadiusArray(triangleStartCounter) = tubeRadiusArray(triangleStartCounter);
        tubeNewRadiusArray(triangleStartCounter+tubeCounter) = tubeRadiusArray(triangleStartCounter+tubeCounter);
        
        if triangleHeightInCells > 0
            % Case 1 - Check the illustration above
            AB = triangleLengthInCells*ds;
            BC = triangleHeightInCells*ds;
            
            for heightCounter = 1:tubeCounter-1
                AD = heightCounter*ds;
                DE = (AD*BC)/AB;
                numCells = round(DE/ds);
                tubeNewRadiusArray(triangleStartCounter+heightCounter)...
                    = tubeRadiusArray(triangleStartCounter+heightCounter)+numCells;
            end
            
        else
            % Case 2 - Check the illustration above
            AB = triangleLengthInCells*ds;
            BC = abs(triangleHeightInCells)*ds;
            
            for heightCounter = 1:tubeCounter-1
                AD = (tubeCounter-heightCounter)*ds;
                DE = (AD*BC)/AB;
                numCells = round(DE/ds);
                
                tubeNewRadiusArray(triangleStartCounter+heightCounter)...
                    = tubeRadiusArray(triangleStartCounter+tubeCounter)+numCells;
            end
        end
        
        % Set the triangleStartCounter
        triangleStartCounter = triangleStartCounter+tubeCounter;
    end
end