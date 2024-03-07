
%% Init

    clear
    close all
    clc
    addpath(genpath('SourceCode'))

%%   Load the surface and prepare the surface normals
    
    nFreqsSim = 14;
    structMeshPreparation = struct();
    structMeshPreparation.orientation = [ 0 0 0];
    structMeshPreparation.position = [ 0 0 0];
    structMeshPreparation.vertexScaling = 1/1000;
    structMeshPreparation.FLIPNORMALS = 1;
    structMeshPreparation.vecFreqSim = linspace( 20e3, 85e3, nFreqsSim ); 
    structMeshPreparation.fileNameMesh = "Data/Models/wall.stl";

    structMeshPreparation.BRDFTransitionPosition = 0.4;
    structMeshPreparation.BRDFTransitionSlope = 2;
    structMeshPreparation.BRDFExponentSpecular = linspace( 8,5,nFreqsSim);
    structMeshPreparation.BRDFExponentDiffractive = linspace( 70,70,nFreqsSim);

    structMeshPreparation.materialStrengthSpecular = 10*linspace( 1,0.8,nFreqsSim);
    structMeshPreparation.materialStrengthDiffractive = 0.05*linspace( 0.5,1,nFreqsSim);
    structMeshPreparation.materialSTransitionPosition = structMeshPreparation.BRDFTransitionPosition;
    structMeshPreparation.materialSTransitionSlope = 2;
    
    structMeshPreparation.precomputeCurvature = 1;
    
    structSurface = prepareMeshSurface( structMeshPreparation, 1 );

%% Setup the structs for processing

    % Make the receiver arrays
    arraySingleEar =  generateCircularArray( 0.02, 0.01, 0);
    nMicsSingleEar = size( arraySingleEar, 1 );
    arrayLeft = [ zeros( nMicsSingleEar, 1 )  arraySingleEar ] + [ -0.010 0.01 0.005];
    arrayRight = [ zeros( nMicsSingleEar, 1 )  arraySingleEar ] + [ -0.010 -0.01 0.005];
    coordsReceivers = [ arrayLeft ; arrayRight ];
    coordsEmitter = [ 0 0 0 ];
    idxArrayLeft = 1 : nMicsSingleEar;
    idxArrayRight = (nMicsSingleEar+1) : (2*nMicsSingleEar);
    
    % Combine it all into the sensor
    structSensor = struct();
    structSensor.position = [ 1 0 0]';
    structSensor.orientation = [ 0 0 180];
    structSensor.coordsEmitter = [ 0 0 -0.01];
    structSensor.coordsReceivers = coordsReceivers;
    structSensor.nMics = size( structSensor.coordsReceivers, 1 );

    % Struct for the parameters of the simulation
    structSimulationParameters = struct();
    structSimulationParameters.doPlot = 0;
    structSimulationParameters.numSamplesImpresp = 16000;
    structSimulationParameters.sampleRateImpresp = 450e3;
    structSimulationParameters.limitsAzimuth = [-20 20];
    structSimulationParameters.limitsElevation = [-20 20];
    structSimulationParameters.numberOfDirections = 50000;
    structSimulationParameters.numberOfDirectionsPerCall = 150000;
    structSimulationParameters.vecFreqSim = structMeshPreparation.vecFreqSim;
    structSimulationParameters.numSamplesIRFilter = 256;
    structSimulationParameters.IRFilterGaussAlpha = 5;
    structSimulationParameters.numDiffractionPoints = 5000;
    structSimulationParameters.approximateImpulseResponseCutDB = -90;
    structSimulationParameters.approximateImpulseResponse = 0;
    structSimulationParameters.ditherRaytracing = 1;
    structSimulationParameters.speedOfSound = 343;

%% Plot the setup, to make sure normals are OK

    nNormals = size( structSurface.surfaceNormals, 1 );
    normalSpacing = round( nNormals / 1000);
    idxNormalsPlot = 1 : normalSpacing : nNormals;
    figure(1212)       
        clf
        hp = patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 1 0 0], 'EdgeAlpha', 0.4); 
        xlabel( 'X-Axis' )
        ylabel( 'Y-Axis' )
        zlabel( 'Z-Axis' )
        grid on
        
        hold on;
            drawTriad( structSensor.position(:), structSensor.orientation(:), 0.05)
        hold off;

        hold on
        for cntNormal = 1 : length( idxNormalsPlot )
            curIdxPlot = idxNormalsPlot( cntNormal );
            curPointOrigin = structSurface.surfaceVertices( curIdxPlot, : );
            curPointNormal = curPointOrigin + structSurface.surfaceNormals( curIdxPlot, : ) * 0.02;
            curPointsPlot = [ curPointOrigin ; curPointNormal ];
            plot3( curPointsPlot(:,1), curPointsPlot(:,2), curPointsPlot(:,3), 'g' );
        end
        hold off;
        axis equal
        set( gca, 'view', [21.6026   19.6988]);

%% Process the raytracing
  
    tic
    structSimulationResult = calculateImpulseResponseFast( structSensor, structSurface, structSimulationParameters );
    toc

%% Process the impulseresponses  

    sigEmit = fm_sweep( 150e3, 40e3, 450e3, 0.4, 1, 10 );
    dataMics = zeros( size( structSimulationResult.impulseResponse ) );
    for cntMic = 1 : size(dataMics,2)
        dataMics( :, cntMic ) = conv( structSimulationResult.impulseResponse(:, cntMic), sigEmit, 'same' );
    end

%% Plot the final result, including the reflection points.

    strengthsReflections = structSimulationResult.strengthsReflected; 
    pointsReflections = structSimulationResult.pointsReflected;
    nPointsReflected = size( pointsReflections, 1 );
    plotPointsSkip = 10;

    idxFrequencyPoints = 1;
    idxMicPoints = 1;
    strengthsReflectionsNorm = normLin( squeeze( strengthsReflections( :, idxFrequencyPoints, idxMicPoints ) ) ); 

    dbCutoffPlot = 60;
    strengthsReflectionsNorm = normLog( strengthsReflectionsNorm, - dbCutoffPlot );
    strengthsReflectionsNorm = normLin( strengthsReflectionsNorm + dbCutoffPlot );
    
    strengthsReflectionsIndexer = round( strengthsReflectionsNorm *254 ) +1;
    cmapPlot = parula(256);
    
    colorsPoints = cmapPlot( strengthsReflectionsIndexer, : );
        figure()       
        
            subplot( 2,3,[1 2 4 5])
                cla
                hp = patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 0.2 0.3 0.3 ], 'EdgeAlpha', 0.3); 
                xlabel( 'X-Axis' )
                ylabel( 'Y-Axis' )
                zlabel( 'Z-Axis' )
                grid on
                set( gca, 'view', [67.1313 7.3250] );
                axis equal

                hold on;
                    drawTriad( structSensor.position(:), structSensor.orientation(:), 0.01)
                    for cntPoint = 1 : plotPointsSkip: nPointsReflected
                        plot3( pointsReflections(cntPoint, 1), pointsReflections(cntPoint,2), pointsReflections(cntPoint,3), '.', 'Color', colorsPoints(cntPoint,:), 'MarkerSize', strengthsReflectionsNorm(cntPoint)*20+0.1 )
                    end
                hold off;

                % Create a colorbar
                c = colorbar;
                c.Label.String = 'Normalized Reflection Strength based on BRDF (dB)';
                colormap(cmapPlot);
                caxis([-dbCutoffPlot, 0]);
                colorbarPosition = [0.1, 0.1, 0.02, 0.8];
                set(c, 'Position', colorbarPosition);
                c.Label.FontSize = 14;  % Adjust the font size as needed
                c.FontSize = 12;
                set( gcf, 'position', [273         348        1725         990])
                title('Reflectivity plot of scene with sensor position/orientation', 'fontsize', 16)

                camlight
                hp.FaceLighting = 'gouraud';
                hp.AmbientStrength = 0.4;
                hp.DiffuseStrength = 0.8;
                hp.SpecularStrength = 0.1;
                hp.SpecularExponent = 1;
                hp.BackFaceLighting = 'unlit';
                
            
            numSamplesPlotIR = 15000;
            IRLeft =  structSimulationResult.impulseResponse(1:numSamplesPlotIR,1);
            IRRight =  structSimulationResult.impulseResponse(1:numSamplesPlotIR,2);
            timeVec = (1:numSamplesPlotIR)/structSimulationParameters.sampleRateImpresp;
            subplot(3,3,3)
                plot( timeVec*1000,IRLeft, 'linewidth', 1.5 )
                grid on
                xlabel( 'Time (ms)')
                ylabel( 'Impulse Response (pressure - au)')
                title( 'Left Impulse Response' )
            subplot(3,3,6)
                plot( timeVec*1000,IRRight, 'linewidth', 1.5 )                
                grid on
                xlabel( 'Time (ms)')
                ylabel( 'Impulse Response (pressure - au)')
                title( 'Right Impulse Response' )
                
            subplot(3,3,9)
                plot( timeVec*1000,IRRight, 'linewidth', 1.5 ) 
                hold on
                plot( timeVec*1000,IRLeft, 'linewidth', 1.5 )
                hold off;
                grid on
                xlabel( 'Time (ms)')
                ylabel( 'Impulse Response (pressure - au)')
                title( 'Impulse Response Overlay' )              
                legend( 'Right', 'Left')

