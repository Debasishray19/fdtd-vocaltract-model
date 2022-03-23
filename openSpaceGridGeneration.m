function [simGridParam, PV_N, listenerInfo] = openSpaceGridGeneration(gridCellTypes, pmlSwitch, pmlLayer)

    % Define domain size
    simGridParam.domainX = 20;
    fprintf("Enter size of the domain along X= %d\n", simGridParam.domainX);
    
    simGridParam.domainY = 20;
    fprintf("Enter size of the domain along Y= %d\n", simGridParam.domainY);
    
    simGridParam.domainZ = 20;
    fprintf("Enter size of the domain along Z= %d\n", simGridParam.domainZ);
    
    % Derive frame size
    simGridParam.frameX = simGridParam.domainX + 2; %2 = dead cells
    simGridParam.frameX = simGridParam.frameX + 2*pmlSwitch*pmlLayer;
    
    simGridParam.frameY = simGridParam.domainY + 2; %2 = dead cells
    simGridParam.frameY = simGridParam.frameY + 2*pmlSwitch*pmlLayer;
    
    simGridParam.frameZ = simGridParam.domainZ + 2; %2 = dead cells
    simGridParam.frameZ = simGridParam.frameZ + 2*pmlSwitch*pmlLayer;
    
    % Define PV_N array to store pressure and velocity components
    % PV_N(:,:,1) = To store cell pressure
    % PV_N(:,:,2) = To store Vx
    % PV_N(:,:,3) = To store Vy
    % PV_N(:,:,4) = To store Vz
    % PV_N(:,:,5) = To store grid cell types
    PV_N = zeros(simGridParam.frameY, simGridParam.frameX, simGridParam.frameZ, 5);
    
    % Define cell types and store it in PV_N(,,5)
    % Declare all the cells as air by default
    PV_N(1:simGridParam.frameY, 1:simGridParam.frameX, 1:simGridParam.frameZ, 5) = gridCellTypes.cell_air;
    
    % Find the mid point of the domain for open-space simulation
    srcX = ceil(simGridParam.frameX/2);
    srcY = ceil(simGridParam.frameY/2);
    srcZ = ceil(simGridParam.frameZ/2);

    % Define the source excitation position to induce energy
    PV_N(srcY, srcX, srcZ, 5) = gridCellTypes.cell_excitation;
    
    wallCondition = 1;
    fprintf('Activate wall inside the  open space [1-Yes 2-No] = %d\n', wallCondition);
    
    if wallCondition == 1
        PV_N(2:simGridParam.frameY-1, srcX+2, 2:simGridParam.frameZ-1, 5) = gridCellTypes.cell_wall;
    end
    
    listenerInfo.listenerX = srcX+5;
    listenerInfo.listenerY = srcY;
    listenerInfo.listenerZ = srcZ;    
end