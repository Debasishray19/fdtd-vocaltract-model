%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE UNITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meter = 1;
centimeter  =1e-2 * meter;

second    = 1;
milisecond = 1e-3 * second;
hertz     = 1/second;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;

gram      = 1e-3;
kilogram  = 1e3*gram;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = 1.140*kilogram/(meter^3);      % Air density  [kg/m^3]
c   = 350*meter/second;              % Sound speed in air [m/s]
maxSigmadt = 0.5;                    % Attenuation coefficient at the PML layer
srate = 44100*srate_mul*hertz;       % Sample frequency
pmlLayer = 6;                        % Number of PML layers
baffleSwitchFlag =(baffleSwitch==1); % By default model should not have head/circular baffle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/srate;                        % Temporal resolution/ sample time period
dx = dt*c*sqrt( 3.0 );               % Spatial resolution along x-direction: CFL Condition
dy = dt*c*sqrt( 3.0 );               % Spatial resolution along x-direction: CFL Condition
dz = dt*c*sqrt( 3.0 );               % Spatial resolution along x-direction: CFL Condition
AudioTime = dur*second;              % Total audio signal time
kappa = rho*c*c;                     % Bulk modulus
ds = dx;                             % Spatial resolution(ds) = dx = dy
rhoC = rho*c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:dt:AudioTime-dt;            % time steps
STEPS = length(t);                % Total time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE CELL TYPES, BETA AND SIGMAPRIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridCellTypes.cell_wall       = 0;
gridCellTypes.cell_air        = 1;
gridCellTypes.cell_excitation = 2;
gridCellTypes.cell_pml0       = 3;
gridCellTypes.cell_pml1       = 4;
gridCellTypes.cell_pml2       = 5;
gridCellTypes.cell_pml3       = 6;
gridCellTypes.cell_pml4       = 7;
gridCellTypes.cell_pml5       = 8;
gridCellTypes.cell_dynamic    = 9;
gridCellTypes.cell_dead       = 10;
gridCellTypes.cell_noPressure = 11;
gridCellTypes.cell_head       = 12;
cell_numTypes   = 13;

vis_Boundary = 2000;
sigmadt = zeros(pmlLayer, 1);

% To store beta(tube wall) and sigmaPrimedt(PML Layers) [beta, sigmaPrimedt]. 
% beta_air = 1, beta_PML = 1 and beta_wall = 0
% sigmaPrimedt = sigmaPrime*dt
% sigmaPrime = 1 - beta + sigma
% e.g - 
% sigma=0 for all the non-PML layers. Hence, sigmaPrime = 1 - beta
% inside the domain. Therefore,
% WALL -> beta = 0, sigma_prima*dt = (1-0)*dt = 1*dt = dt 
% AIR  -> beta = 1, sigma_prima*dt = (1-1)*dt = 0*dt = 0
% HEAD CELLS -> beta = 0, sigma_prima*dt = (1-0)*dt = 1*dt = dt 
% [NOTE] - We are considering excitation cells/head cells as special wall cells
typeValues = zeros(2, cell_numTypes);
typeValues(:, gridCellTypes.cell_wall+1) = [0, dt];       % VT walls
typeValues(:, gridCellTypes.cell_air+1) = [1, 0];         % air
typeValues(:, gridCellTypes.cell_noPressure+1) = [1, 0];  % air
typeValues(:, gridCellTypes.cell_excitation+1) = [0, dt]; % excitation
typeValues(:, gridCellTypes.cell_head+1) = [0, dt];       % head cells

% Define beta and sigmaPrimedt for PML layers
% For PML layers beta = 1,
% sigmaPrimedt = [1-beta+sigma]*dt = [1-1+sigma]*dt = sigmadt
for pmlCounter = 0:pmlLayer-1
    sigmadt(pmlCounter+1) = (pmlCounter/(pmlLayer-1)) * maxSigmadt;
    typeValues(:, gridCellTypes.cell_pml0+1+pmlCounter) = [1, sigmadt(pmlCounter+1)];
end
typeValues(:, gridCellTypes.cell_dead+1) = [0, 1000000];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION TYPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Note]: I have created two separate functions for the circular and 
% elliptical tube generation. However, that's not required since circle is
% a special case of ellipse where the ratio between the semiMajorAxis and
% semiMinorAxis is equal. [semiMajorAxis:semiMinorAxis = 1:1]

switch simulationType
    case 1 %Simulation: Wave propagation from a point source inside an open space       
        [simGridParam, PV_N, listenerInfo] = openSpaceGridGeneration(gridCellTypes, pmlSwitch, pmlLayer);
        
    case 2 %Simulation: Wave propagation inside a regular cylindrical tube
        if cross_sectionType == 1
            [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo]...
                = vt_circularTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds);
        elseif cross_sectionType == 2
            [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo]...
                = vt_ellipticalTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds);
        elseif cross_sectionType == 3
            [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo]...
                = vt_squareTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds);
        end
                   
    case 3 %Simulation: Wave propagation inside a tube having VT geometry
        if cross_sectionType == 1
            [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo]...
                = vt_circularTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds);
        elseif cross_sectionType == 2
            [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo]...
                = vt_ellipticalTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds);
        elseif cross_sectionType == 3
            [simGridParam, PV_N, tubeStartArea, tubeStart, tubeEnd, totalTubeLengthInCells, currTubeSectionDiameterCells_SegmentCounter, listenerInfo, sourceInfo]...
                = vt_squareTubeGeneration(simulationType, junctionType, vowel, boundaryInterpolation, gridCellTypes, baffleSwitchFlag, pmlSwitch, pmlLayer, rad, ds);
        end
        
    otherwise
        fprintf("Incorrect input. Rerun the simulation");
        return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VARIABLES TO SAVE PRESSURE AT MOUTH_END AND GLOTTAL_END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pr_mouthEnd = zeros(1,STEPS);
Pr_glottalEnd = zeros(1,STEPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETERMINE WALL IMPEDANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu3D = 0.005; %Boundary Admittance Coefficient
alpha3D = 1/(0.5+0.25*(mu3D +(1/mu3D))); % Sound Absorption Coefficient
z_inv = 1 / (rhoC*((1+sqrt(1-alpha3D))/(1-sqrt(1-alpha3D))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOURCE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sourceModelType = sourceType; 

switch sourceModelType
    case 1 % Sine wave source model
        excitationF = 440*hertz;            
        srcAmplitude = 25;
        exeT = linspace(1, STEPS, STEPS);
        excitationV = srcAmplitude * sin(2*pi*excitationF*dt*(exeT(:)-1));

    case 2 % Gaussian source model
        f0 = 10*kilohertz;
        bellPeakPos = 0.646/f0;
        bellWidth = 0.29*bellPeakPos;
        excitationV = exp(-((t-bellPeakPos)./bellWidth).^2);
        
    case 3 % Impulse response function
        excitationV = src_ImpulseSignal(srate, 10000, 2, 22050);
        
    case 4 % gaussian white noise
        excitationV = randn(STEPS, 1); %VIC this will be regenerated in each simulation
        
    case 5 % Vocal fold model - Two Mass Model
        % Set the vocal fold parameters
        [airParam, vf_structuralParam, vf_flowParam, vf_matParam] = vf_SetVocalFoldParams();
    otherwise
        fprintf("Incorrect input. Rerun the simulation");
        return;     
end

% Define source propagation direction
% srcDirection index mean: 1 = Left  = -1
%                          2 = Down  = -1
%                          3 = Back  = -1
%                          4 = Right =  1
%                          5 = Up    =  1
%                          6 = Front =  1
srcDirection = [0 0 0 1 0 0]; % For forward direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOURCE AND VIRTUAL MIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frameX = simGridParam.frameX;
frameY = simGridParam.frameY;
frameZ = simGridParam.frameZ;

% Find the mid point of the domain for open-space simulation
midX = floor(frameX/2);
midY = floor(frameY/2);
midZ = floor(frameZ/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE PML LAYERS AND DEAD CELLS INSIDE THE 3D GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define cell_dead to the outer most layer of the frame

%---Adding dead cells to the top and bottom surfaces---%
PV_N(1, 1:frameX, 1:frameZ,5) = gridCellTypes.cell_dead;
PV_N(frameY, 1:frameX, 1:frameZ,5) = gridCellTypes.cell_dead;

%---Adding dead cells to the back and front surfaces---%
PV_N(1:frameY, 1:frameX, 1,5) = gridCellTypes.cell_dead;
PV_N(1:frameY, 1:frameX, frameZ,5) = gridCellTypes.cell_dead;

%---Adding dead cells to the left and right surfaces---%
PV_N(1:frameY, 1, 1:frameZ,5) = gridCellTypes.cell_dead;
PV_N(1:frameY, frameX, 1:frameZ,5) = gridCellTypes.cell_dead;

if pmlSwitch == 1
  
   % Define PML layers starting from the outermost layer
   pmlType = gridCellTypes.cell_pml5;
   
   xShift = 1;
   xStart = 2;
   xEnd   = frameX-1;
   
   yShift = 1;
   yStart = 2;
   yEnd   = frameY-1;
   
   zShift = 1;
   zStart = 2;
   zEnd   = frameZ-1;
   
   for pmlCount=1:pmlLayer
       %---Adding pml layers to the top and bottom surfaces---%
       PV_N(yShift+pmlCount, xStart:xEnd, zStart:zEnd,5) = pmlType;
       PV_N(frameY-pmlCount, xStart:xEnd, zStart:zEnd,5) = pmlType;
       
       %---Adding pml layers to the back and front surfaces---%
       PV_N(yStart:yEnd, xStart:xEnd, zShift+pmlCount,5) = pmlType;
       PV_N(yStart:yEnd, xStart:xEnd, frameZ-pmlCount,5) = pmlType;
       
       %---Adding pml layers to the left and right surfaces---%
       PV_N(yStart:yEnd, xShift+pmlCount, zStart:zEnd,5) = pmlType;
       PV_N(yStart:yEnd, frameX-pmlCount, zStart:zEnd,5) = pmlType;
       
       % Update the pmlLayer type.
       % Update the start and end position
       pmlType = pmlType-1;
       xStart = xStart+1;
       xEnd = xEnd-1;
       
       yStart = yStart+1;
       yEnd = yEnd-1;
       
       zStart = zStart+1;
       zEnd = zEnd-1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL VISUALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualize the PML layers and excitation cell by taking a slice of the domain
buildFrame = PV_N(:,:,:,5);
%buildFrame(midY, midX, midZ) = -1; % Uncomment it to visualize source position for the openspace simulation 
buildFrame(listenerInfo.listenerY, listenerInfo.listenerX, listenerInfo.listenerZ) = -1;

% Slice along XZ-plane
frame = buildFrame(ceil(frameY/2),:,:);
frame_vis1_xz = reshape(frame, [frameX,frameZ]);
figure;
imagesc(frame_vis1_xz);

% Slice along YZ-plane
frame = buildFrame(:,ceil(frameX/2),:);
frame_vis2_yz = reshape(frame, [frameY,frameZ]);
figure;
imagesc(frame_vis2_yz);

% Slice along XY-plane
frame = buildFrame(:,:,ceil(frameZ/2));
frame_vis3_xy = reshape(frame, [frameY,frameX]);
figure;
imagesc(frame_vis3_xy);

wavePropagationVis = zeros(frameY, frameX, frameZ);
vis_Boundary = 2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOK UP TABLE FOR BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wall = gridCellTypes.cell_wall;
Air = gridCellTypes.cell_air;

% Look up table store cell types in the following order: 
% [right_current, up_current, front_current]
lookUp_table = [Air, Air, Air;...
                Air, Air, Wall;...
                Air, Wall, Air;...
                Air, Wall, Wall;...
                Wall, Air, Air;...
                Wall, Air, Wall;...
                Wall, Wall, Air;...
                Wall, Wall, Wall];

air_normalV_component = [0, 1, 1, 0.7071, 1, 0.7071, 0.7071, 0.5774];
wall_normalV_component = [0.5774, 1, 1, 0.7071, 1, 0.7071, 0.7071, 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MINBETA AND MAXSIGMAPRIMEDT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -> Store the cellTypes, typeIndex, beta and sigmaPrimedt separately
% -> Compute minbeta and maxsigmaprimedt

height = frameY-2; % => Number of cells along Y - cell_Dead
width  = frameX-2; % => Number of cells along X - cell_Dead
depth  = frameZ-2; % => Number of cells along Z - cell_Dead

% Create arrays to store cellTypes, beta and sigmaPrime_dt values
% We create 4 columns here as each grid cell has three neighbors
% Note=> We don't need to define these arrays [Just for my understanding]
cellTypes     = zeros(1, 4);
typeIndex     = zeros(1, 4);
beta          = zeros(1, 4);
sigmaPrime_dt = zeros(1, 4);

% To store is_excitation and excitation_weight
numGridCellsForComputation = (frameX-2)*(frameY-2)*(frameZ-2);
is_excitation              = zeros(4, numGridCellsForComputation);
are_we_not_excitations     = zeros(1, 3);
excitation_weight          = zeros(3, numGridCellsForComputation);
xor_term                   = zeros(6, numGridCellsForComputation);
N_out                      = zeros(1, numGridCellsForComputation);
N_in                       = zeros(1, numGridCellsForComputation);
% Create arrays to store minBeta and maxSigmaPrime_dt
% Note=> We consider minBeta, because beta parameter is defined at the centre 
% of the ccell. And the velocity components are defined either on edges [2D
% grid cell]/ side surfaces [3D grid cell]. Therefore, we tie-break between
% these grid cells considering minBeta value.

minVxBeta = zeros(height, width, depth);
minVyBeta = zeros(height, width, depth);
minVzBeta = zeros(height, width, depth);

maxVxSigmaPrimedt = zeros(height, width, depth);
maxVySigmaPrimedt = zeros(height, width, depth);
maxVzSigmaPrimedt = zeros(height, width, depth);

% Store the excitation_weight and are_we_not_excitations in matrix format 
% to compute Vx, Vy, Vz
excitation_Vx_weight = zeros(height, width, depth);
excitation_Vy_weight = zeros(height, width, depth);
excitation_Vz_weight = zeros(height, width, depth);

are_we_not_excitations_Vx = zeros(height, width, depth);
are_we_not_excitations_Vy = zeros(height, width, depth);
are_we_not_excitations_Vz = zeros(height, width, depth);

% Store the xor_val in matrix format
xor_val1 = zeros(height, width, depth);
xor_val2 = zeros(height, width, depth);
xor_val3 = zeros(height, width, depth);
xor_val4 = zeros(height, width, depth);
xor_val5 = zeros(height, width, depth);
xor_val6 = zeros(height, width, depth);

% Store the N_out and N_in in matrix format
N_out_mat = zeros(height, width, depth);
N_in_mat  = zeros(height, width, depth);

% Note=> We are reshaping sigmaPrimedt in matrix format as we need this
% while calculating pressure [Check the denominator of the discretized pressure equation]
sigmaPrimedt = zeros(height, width, depth);

% Note=> Here for each grid cell, we store cell type, beta and siggmaPrine_dt
% values in arrays for its neighbor cells and the current cell. In case of
% For example, in case of - 
% 2D/2.5D => [current, right_current, top_current]
% 3D => [current, right_current, up_current, front_current]

% Set a counter to save data for is_excitation, excitation_weight, xor,
% N_out, N_in
counter = 1;

for height_idx = 2:frameY-1
    for width_idx = 2:frameX-1
        for depth_idx = 2:frameZ-1
            
            % Find the cellTypes = [current, right_current, up_current, front_current]
            cellTypes = ...
                [PV_N(height_idx, width_idx, depth_idx, 5), PV_N(height_idx, width_idx+1, depth_idx, 5), PV_N(height_idx-1, width_idx, depth_idx, 5), PV_N(height_idx, width_idx, depth_idx+1, 5)];
            
            % For typeIndex add 1 to cellTypes
            typeIndex = cellTypes + 1;
            
            % Store beta values in beta array
            beta = typeValues(1, typeIndex(:)); 
            
            % Calculate minBeta
            min_beta_Vx = min(beta([1,2])); % minBeta(current, right_current)
            min_beta_Vy = min(beta([1,3])); % minBeta(current, top_current)
            min_beta_Vz = min(beta([1,4])); % minBeta(current, front_current
            
            minVxBeta(height_idx-1, width_idx-1, depth_idx-1) = min_beta_Vx;
            minVyBeta(height_idx-1, width_idx-1, depth_idx-1) = min_beta_Vy;
            minVzBeta(height_idx-1, width_idx-1, depth_idx-1) = min_beta_Vz;
            
            % Store sigmaPrime*dt values in sigmaPrime_dt
            sigmaPrime_dt = typeValues(2, typeIndex(:));
                        
            % Store sigmaPrime_dt of the current cell only
            sigmaPrimedt(height_idx-1, width_idx-1, depth_idx-1)  = sigmaPrime_dt(1);
            
            % Calculate maxSigmaPrime_dt
            max_sigmaPrimedt_Vx = max(sigmaPrime_dt([1,2])); % maxsigmaPrimedt(current, right_current)
            max_sigmaPrimedt_Vy = max(sigmaPrime_dt([1,3])); % maxsigmaPrimedt(current, top_current)
            max_sigmaPrimedt_Vz = max(sigmaPrime_dt([1,4])); % maxsigmaPrimedt(current, front_current)
            
            maxVxSigmaPrimedt(height_idx-1, width_idx-1, depth_idx-1) = max_sigmaPrimedt_Vx;
            maxVySigmaPrimedt(height_idx-1, width_idx-1, depth_idx-1) = max_sigmaPrimedt_Vy;
            maxVzSigmaPrimedt(height_idx-1, width_idx-1, depth_idx-1) = max_sigmaPrimedt_Vz;
            
            % Check whether the current cell is an excitation cell or not
            is_excitation(:, counter) = [cellTypes(1) == gridCellTypes.cell_excitation, cellTypes(2) == gridCellTypes.cell_excitation, cellTypes(3) == gridCellTypes.cell_excitation, cellTypes(4) == gridCellTypes.cell_excitation];
            
            % Verify if both the current cell and neighboring cells are not excitation cells
            are_we_not_excitations = [(1-is_excitation(1, counter)) .* (1-is_excitation(2, counter)),...
                                      (1-is_excitation(1, counter)) .* (1-is_excitation(3, counter)),...
                                      (1-is_excitation(1, counter)) .* (1-is_excitation(4, counter))];
            
            are_we_not_excitations_Vx (height_idx-1, width_idx-1, depth_idx-1) = are_we_not_excitations(1);
            are_we_not_excitations_Vy (height_idx-1, width_idx-1, depth_idx-1) = are_we_not_excitations(2); 
            are_we_not_excitations_Vz (height_idx-1, width_idx-1, depth_idx-1) = are_we_not_excitations(3); 
            
            % Determine excitation weight for the current cell
            % Find excitation velocity going out of the cell and coming
            % back into the cells depending upon the source direction. 
            % Then add them all to find the net velocity [velocity components]
            % for the current cell.
               
            excitation_weight_forward = is_excitation(1, counter).*(srcDirection(4:6)');
            excitation_weight_backward = is_excitation(2:4, counter).*(srcDirection(1:3)');
            excitation_weight(:, counter) = excitation_weight_forward + excitation_weight_backward;
            
            % Store the excitation_weight in matrix format
            excitation_Vx_weight(height_idx-1, width_idx-1, depth_idx-1) = excitation_weight(1, counter);
            excitation_Vy_weight(height_idx-1, width_idx-1, depth_idx-1) = excitation_weight(2, counter);
            excitation_Vz_weight(height_idx-1, width_idx-1, depth_idx-1) = excitation_weight(3, counter);
            
            % Check the adjacent cells to determine the velocity
            % direction
            xor_val = [beta(2) .* (1-beta(1)), beta(1) .* (1-beta(2)), ...
                       beta(3) .* (1-beta(1)), beta(1) .* (1-beta(3)),...
                       beta(4) .* (1-beta(1)), beta(1) .* (1-beta(4))];
                   
            xor_val1(height_idx-1, width_idx-1, depth_idx-1) = xor_val(1);
            xor_val2(height_idx-1, width_idx-1, depth_idx-1) = xor_val(2);
            xor_val3(height_idx-1, width_idx-1, depth_idx-1) = xor_val(3);
            xor_val4(height_idx-1, width_idx-1, depth_idx-1) = xor_val(4);
            xor_val5(height_idx-1, width_idx-1, depth_idx-1) = xor_val(5);
            xor_val6(height_idx-1, width_idx-1, depth_idx-1) = xor_val(6);
            
            % Find the boundary condition from the lookup table and
            % select the correspoding unit vector
            [isbetaMatched,lookup_idx] = ismember(beta(2:4), lookUp_table, 'rows');
               
            N_out(1, counter) = air_normalV_component(lookup_idx)*beta(1);
            N_in(1, counter)  = wall_normalV_component(lookup_idx)*(1-beta(1));
            
            N_out_mat(height_idx-1, width_idx-1, depth_idx-1) = N_out(1, counter);
            N_in_mat(height_idx-1, width_idx-1, depth_idx-1)  = N_in(1, counter);
            
            % Increment the counter
            counter = counter+1;
        end
    end
end

betaVxSqr = minVxBeta.*minVxBeta;
betaVxSqr_dt_invRho = (betaVxSqr*dt)/rho;

betaVySqr = minVyBeta.*minVyBeta;
betaVySqr_dt_invRho = (betaVySqr*dt)/rho;

betaVzSqr = minVzBeta.*minVzBeta;
betaVzSqr_dt_invRho = (betaVzSqr*dt)/rho;

rho_sqrC_dt_invds = (kappa*dt)/ds;
rho_sqrC_dt       = kappa*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISE ACOUSTIC PARAMETERS AND COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV_NPlus1  =  zeros(frameY, frameX, frameZ, 5);

CxVx = zeros(frameY-2, frameX-2, frameZ-2);
CyVy = zeros(frameY-2, frameX-2, frameZ-2);
CzVz = zeros(frameY-2, frameX-2, frameZ-2);

CxP  = zeros(frameY-2, frameX-2, frameZ-2);
CyP  = zeros(frameY-2, frameX-2, frameZ-2);
CzP  = zeros(frameY-2, frameX-2, frameZ-2);

Pr_next = zeros(frameY-2, frameX-2, frameZ-2);
Vx_next = zeros(frameY-2, frameX-2, frameZ-2);
Vy_next = zeros(frameY-2, frameX-2, frameZ-2);
Vz_next = zeros(frameY-2, frameX-2, frameZ-2);

vb_alphaX = zeros(frameY-2, frameX-2, frameZ-2);
vb_alphaY = zeros(frameY-2, frameX-2, frameZ-2);
vb_alphaZ = zeros(frameY-2, frameX-2, frameZ-2);

% STEPXX: Introduce Dirichlet Boundary Condition
if rad==2 && simulationType~=1
    xNoPressure = tubeStart.startX + totalTubeLengthInCells;
    gridNoPressurePlane = PV_N(:, xNoPressure, :, 5);
    isNoPressureCell = (gridNoPressurePlane~=gridCellTypes.cell_noPressure);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD Acceleration using Parallel Computing Toolbox and GPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCTimplementation = 0;

if PCTimplementation == 1
    CxVx = gpuArray(CxVx);
    CyVy = gpuArray(CyVy);
    CzVz = gpuArray(CzVz);
    CxP = gpuArray(CxP);
    CyP = gpuArray(CyP);
    CzP = gpuArray(CzP);
    
    Pr_next = gpuArray(Pr_next);
    Vx_next = gpuArray(Vx_next);
    Vy_next = gpuArray(Vy_next);
    Vz_next = gpuArray(Vz_next);
    PV_N = gpuArray(PV_N);
    PV_NPlus1 = gpuArray(PV_NPlus1);
    
    isNoPressureCell = gpuArray(isNoPressureCell);
    excitationV = gpuArray(excitationV);
    
    rho_sqrC_dt_invds = gpuArray(rho_sqrC_dt_invds);
    minVxBeta = gpuArray(minVxBeta);
    minVyBeta = gpuArray(minVyBeta);
    minVzBeta = gpuArray(minVzBeta);
    maxVxSigmaPrimedt = gpuArray(maxVxSigmaPrimedt);
    maxVySigmaPrimedt = gpuArray(maxVySigmaPrimedt);
    maxVzSigmaPrimedt = gpuArray(maxVzSigmaPrimedt);
    
    betaVxSqr_dt_invRho = gpuArray(betaVxSqr_dt_invRho);
    betaVySqr_dt_invRho = gpuArray(betaVySqr_dt_invRho);
    betaVzSqr_dt_invRho = gpuArray(betaVzSqr_dt_invRho);
    excitation_Vx_weight = gpuArray(excitation_Vx_weight);
    excitation_Vy_weight = gpuArray(excitation_Vy_weight);
    excitation_Vz_weight = gpuArray(excitation_Vz_weight);
    
    are_we_not_excitations_Vx = gpuArray(are_we_not_excitations_Vx);
    are_we_not_excitations_Vy = gpuArray(are_we_not_excitations_Vy);
    are_we_not_excitations_Vz = gpuArray(are_we_not_excitations_Vz);
    
    vb_alphaX = gpuArray(vb_alphaX);
    vb_alphaY = gpuArray(vb_alphaY);
    vb_alphaZ = gpuArray(vb_alphaZ);
    
    xor_val1 = gpuArray(xor_val1);
    xor_val2 = gpuArray(xor_val2);
    xor_val3 = gpuArray(xor_val3);
    xor_val4 = gpuArray(xor_val4);
    xor_val5 = gpuArray(xor_val5);
    xor_val6 = gpuArray(xor_val6);
    
    N_out_mat = gpuArray(N_out_mat);
    N_in_mat = gpuArray(N_in_mat);
    
    z_inv = gpuArray(z_inv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deb: Just to verify if the geometry is correct or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometryArray = [approx. cross-sectional area, , , tube length]

if simulationType~=1
    geometryArray = zeros(4,totalTubeLengthInCells);
    geometryArrayCounter = 1;

    for tubeAreaCounter = tubeStart.startX:tubeStart.startX + totalTubeLengthInCells-1
        [gridPlaneProp, gridCellTypeInplane] = vt_findCellTypes(PV_N, gridCellTypes, tubeAreaCounter);
        cellCountInVTContour = nnz(gridPlaneProp == gridCellTypeInplane.inVTContour);
        geometryArray(1, geometryArrayCounter) = cellCountInVTContour*(dy*dz)*(100*100);
        geometryArray(2, geometryArrayCounter) = currTubeSectionDiameterCells_SegmentCounter(1,geometryArrayCounter);
        geometryArray(3, geometryArrayCounter) = currTubeSectionDiameterCells_SegmentCounter(2,geometryArrayCounter);
        geometryArray(4, geometryArrayCounter) = geometryArrayCounter*dx;

        geometryArrayCounter = geometryArrayCounter+1;  
    end
end

% Open a new figure window to visualize the simulation
figure;

for T = 1:STEPS
    
%     % Start the clock
%     tic;
    
    % STEP1: Calculate [del.V] = [dVx/dx + dVy/dy + dVz/dz]
    % CxVx = dVx, where Vx = velocity along x direction = Vx_curr - Vx_left
    % CyVy = dVy, where Vy = velocity along y direction = Vy_Curr - Vy_down
    % CzVz = dVz, where Vz = velocity along z direction = Vz_curr - Vz_back
    
    CxVx(1:frameY-2, 1:frameX-2, 1:frameZ-2) = PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) - ...
                                               PV_N(2:frameY-1, 1:frameX-2, 2:frameZ-1, 2);
                                   
    CyVy(1:frameY-2, 1:frameX-2, 1:frameZ-2) = PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) - ...
                                               PV_N(3:frameY, 2:frameX-1, 2:frameZ-1, 3);
                                           
    CzVz(1:frameY-2, 1:frameX-2, 1:frameZ-2) = PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) - ...
                                               PV_N(2:frameY-1, 2:frameX-1, 1:frameZ-2, 4);                                          
    
    % STEP2: Calculate Pr_next
    Pr_next(1:frameY-2, 1:frameX-2, 1:frameZ-2) = (PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 1) - (rho_sqrC_dt_invds.*(CxVx+CyVy+CzVz)))./...
        (1+sigmaPrimedt);
    
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 1) = Pr_next(:,:,:);
    
    % Make pressure values of no_pressure cells as zeros
    if simulationType~=1
        PV_NPlus1(:, xNoPressure, :, 1) = PV_NPlus1(:, xNoPressure, :, 1).*isNoPressureCell; 
    end
    
    % STEP3: Calculate Vx_next, Vy_next and Vz_next
    % CxP [del.Px] = [dPx/dx] = Pr_right - Pr_curr
    % CyP [del.Py] = [dPy/dy] = Pr_top   - Pr_curr
    % CzP [del.Pz] = [dPz/dz] = Pr_front - Pr_Curr
    
    CxP = (PV_NPlus1(2:frameY-1, 3:frameX, 2:frameZ-1, 1) - PV_NPlus1 (2:frameY-1, 2:frameX-1, 2:frameZ-1,1))/dx;
    CyP = (PV_NPlus1(1:frameY-2, 2:frameX-1, 2:frameZ-1, 1) - PV_NPlus1 (2:frameY-1, 2:frameX-1, 2:frameZ-1,1))/dy;
    CzP = (PV_NPlus1(2:frameY-1, 2:frameX-1, 3:frameZ, 1) - PV_NPlus1 (2:frameY-1, 2:frameX-1, 2:frameZ-1,1))/dz;
     
    Vx_next(1:frameY-2, 1:frameX-2, 1:frameZ-2) = (minVxBeta.*PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2)) - (betaVxSqr_dt_invRho.*CxP);
    Vy_next(1:frameY-2, 1:frameX-2, 1:frameZ-2) = (minVyBeta.*PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3)) - (betaVySqr_dt_invRho.*CyP);
    Vz_next(1:frameY-2, 1:frameX-2, 1:frameZ-2)	= (minVzBeta.*PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4)) - (betaVzSqr_dt_invRho.*CzP);
    
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) = Vx_next(:,:,:);
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) = Vy_next(:,:,:);
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) = Vz_next(:,:,:);
          
    % STEP4(i) : Inject excitation velocity
    % STEP4(ii): Enforce boundary condition 
    
    exeCurrentVal = excitationV(T);
    
    % Update Vx, Vy and Vz components of the current cell with the
    % excitation velocity
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) = PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) + ...
        exeCurrentVal.*excitation_Vx_weight.*maxVxSigmaPrimedt;
    
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) = PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) + ...
        exeCurrentVal.*excitation_Vy_weight.*maxVySigmaPrimedt;
    
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) = PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) + ...
        exeCurrentVal.*excitation_Vz_weight.*maxVzSigmaPrimedt;
    
    % Determine velocity near the wall
    vb_alphaX = xor_val2.*PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 1).*N_out_mat - ...
                xor_val1.*PV_NPlus1(2:frameY-1, 3:frameX,   2:frameZ-1, 1).*N_in_mat;
            
    vb_alphaY = xor_val4.*PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 1).*N_out_mat - ...
                xor_val3.*PV_NPlus1(1:frameY-2, 2:frameX-1, 2:frameZ-1, 1).*N_in_mat;
            
    vb_alphaZ = xor_val6.*PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 1).*N_out_mat - ...
                xor_val5.*PV_NPlus1(2:frameY-1, 2:frameX-1, 3:frameZ,   1).*N_in_mat;
            
            
    % Apply tube boundary condition/wall losses
    vb_alphaX = (vb_alphaX.*are_we_not_excitations_Vx).*z_inv;
    vb_alphaY = (vb_alphaY.*are_we_not_excitations_Vy).*z_inv;
    vb_alphaZ = (vb_alphaZ.*are_we_not_excitations_Vz).*z_inv;
    
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) = ...
        PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) + maxVxSigmaPrimedt.*vb_alphaX;
                                                   
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) = ...
        PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) + maxVySigmaPrimedt.*vb_alphaY;
                                                   
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) = ...
        PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) + maxVzSigmaPrimedt.*vb_alphaZ;
    
    % Update velocity components
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2) = PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 2)./(minVxBeta+maxVxSigmaPrimedt);
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3) = PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 3)./(minVyBeta+maxVySigmaPrimedt);
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4) = PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 4)./(minVzBeta+maxVzSigmaPrimedt);
    
    % STEP5: Re-store the grid cell type
    PV_NPlus1(2:frameY-1, 2:frameX-1, 2:frameZ-1, 5) = PV_N(2:frameY-1, 2:frameX-1, 2:frameZ-1, 5);
    
    % STEP6: Pass over outer dead cells
    PV_NPlus1(1, :, :, 1:4) = 0;
    PV_NPlus1(1, :, :, 5) = PV_N(1, :, :, 5);
    
    PV_NPlus1(frameY, :, :, 1:4) = 0;
    PV_NPlus1(frameY, :, :, 5) = PV_N(1, :, :, 5);
     
    PV_NPlus1(:, 1, :, 1:4) = 0;
    PV_NPlus1(:, 1, :, 5) = PV_N(:, 1, :, 5);
    
    PV_NPlus1(:, frameX, :, 1:4) = 0;
    PV_NPlus1(:, frameX, :, 5) = PV_N(:, frameX, :, 5); 

    PV_NPlus1(:, :, 1, 1:4) = 0;
    PV_NPlus1(:, :, 1, 5) = PV_N(:, :, 1, 5);
    
    PV_NPlus1(:, :, 1, 1:4) = 0;
    PV_NPlus1(:, :, frameZ, 5) = PV_N(:, :, frameZ, 5); 
    
    % STEP6: Copy PV_Nplus1 to PV_N for the next time step
    PV_N = PV_NPlus1;
    
    % Deb: Print remaining step numbers
    if mod(T,1000)==0
        fprintf('Remaining STEPS = %d\n', STEPS-T);
    end
        
    % Save audio data as change in pressure
    if simulationType~=1
        Pr_mouthEnd(T) = PV_NPlus1(listenerInfo.listenerY, listenerInfo.listenerX, listenerInfo.listenerZ, 1);
        Pr_glottalEnd(T) = PV_NPlus1(sourceInfo.sourceY, sourceInfo.sourceX, sourceInfo.sourceZ, 1);
    end
%     % Stop the clock and record the time
%     stopClock = toc;
%     fprintf('Time taken for the step = %d\n', stopClock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deb: Uncomment the below section to verify the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Take slices of the 3D domain to visualize the wave propagation
    % Slice along xz-plane
%     test_xz= PV_N(ceil(frameY/2),:,:,1);
%     test_2D_xz = reshape(test_xz, [frameX,frameZ]);
%     chech_xz_nan = isnan(test_2D_xz);
    
    % Slice along YZ-plane
%     test_yz = PV_N(:,ceil(frameX/2),:,1);
%     test_2D_yz = reshape(test_yz, [frameY,frameZ]);
%     chech_yz_nan = isnan(test_2D_yz);
    
    % Slice along XY-plane
    test_xy = PV_N(:, :, ceil(frameZ/2),1); % PV_N(:,:,tubeStart.startZ,1);
    test_2D_xy = reshape(test_xy, [frameY,frameX]);
    chech_xy_nan = isnan(test_2D_xy);
 
    wavePropagationVis = test_2D_xy;
    wavePropagationVis(PV_N(:,:,tubeStart.startZ,5) == gridCellTypes.cell_wall) = vis_Boundary; %[comment out this line for open space simulation]
   
    % Plot the pressure wave pripagation
    if plotting==1 && ~mod(T,1)
       imagesc(wavePropagationVis,[-1000 4000]);
       title(['STEP NUMBER: ' num2str(T) ' OUT OF ' num2str(STEPS)]);
       drawnow; 
    end
    
%     if max(max(chech_xz_nan))==1 || max(max(chech_xz_nan))==1 || max(max(chech_xz_nan))==1
%         fprintf('Solver exploded at step = %d', T);
%         return;
%     end  
end

if saveAudioData == 1
    save('audioData.mat', 'Pr_mouthEnd', 'Pr_glottalEnd', 'excitationV', 'srate', 'srate_mul');
end