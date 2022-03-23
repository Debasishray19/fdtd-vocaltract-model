% Generate geometry data
vtModel_geometry = PV_N(:,:,:,5);
simulationDomianSize = size(vtModel_geometry);

nx = simulationDomianSize(2);
ny = simulationDomianSize(1);
nz = simulationDomianSize(3);

% Create the file identifier
geometryFileID = fopen('geometryInfo.txt','w');
fprintf(geometryFileID, '%d\n', nx);
fprintf(geometryFileID, '%d\n', ny);
fprintf(geometryFileID, '%d\n', nz);
count = 0;
for xCount = 1:nx
    for yCount = 1:ny
        for zCount = 1:nz
           count = count+1;
           gridCellTypeValue = vtModel_geometry(yCount, xCount, zCount);
           fprintf(geometryFileID, '%d\n', gridCellTypeValue);
        end  
    end
end

% Close the file identifier
fclose(geometryFileID);

