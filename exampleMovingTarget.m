
%% Init

    clear
    close all
    clc
    addpath(genpath('../../SourceCode'))

%%   Load the surface and prepare the surface normals
    
    nFreqsSim = 16;
    structMeshPreparation = struct();
    structMeshPreparation.orientation = [ 0 0 0];
    structMeshPreparation.position = [ 0 0 0];
    structMeshPreparation.vertexScaling = 1/1000 * 1.5;
    structMeshPreparation.FLIPNORMALS = 1;
    structMeshPreparation.vecFreqSim = linspace( 20e3, 200e3, nFreqsSim ); 
    structMeshPreparation.fileNameMesh = "Data/Models/sphere.stl";

    structMeshPreparation.BRDFTransitionPosition = 0.2;
    structMeshPreparation.BRDFTransitionSlope = 1;
    structMeshPreparation.BRDFExponentSpecular = linspace( 12,4,nFreqsSim);
    structMeshPreparation.BRDFExponentDiffractive = linspace( 70,70,nFreqsSim);
    structMeshPreparation.materialStrengthSpecular = 10*linspace( 1,0.7,nFreqsSim);
    structMeshPreparation.materialStrengthDiffractive = 0.4*linspace( 0.8,1,nFreqsSim);
    structMeshPreparation.materialSTransitionPosition = 0.2;
    structMeshPreparation.materialSTransitionSlope = 1;
    structMeshPreparation.precomputeCurvature = 1;
    
   structSurface = prepareMeshSurface( structMeshPreparation, 1 );

%% Setup the structs for processing

    arraySingleEar =  generateCircularArray( 0.0015, 0.005, 1);
    nMicsSingleEar = size( arraySingleEar, 1 );
    arrayFinal = [ zeros( nMicsSingleEar, 1 )  arraySingleEar ] + [ -0.010 0.01 0.005];
    coordsReceivers = arrayFinal;
    
    structSensor = struct();
    structSensor.position = [ 0.08 0 0];
    structSensor.orientation = [ 0 0 180]';
    structSensor.coordsEmitter = [ 0 0 -0.01];
    structSensor.coordsReceivers = coordsReceivers;
    structSensor.nMics = size( structSensor.coordsReceivers, 1 );
 
    % Struct for the parameters of the simulation
    structSimulationParameters = struct();
    structSimulationParameters.doPlot = 0;
    structSimulationParameters.numSamplesImpresp = 7500;
    structSimulationParameters.sampleRateImpresp = 1e6;
    structSimulationParameters.limitsAzimuth = [-15 15];
    structSimulationParameters.limitsElevation = [-15 15];
    structSimulationParameters.numberOfDirections = 50000;
    structSimulationParameters.numberOfDirectionsPerCall = 150000;
    structSimulationParameters.vecFreqSim = structMeshPreparation.vecFreqSim;
    structSimulationParameters.numSamplesIRFilter = 256;
    structSimulationParameters.IRFilterGaussAlpha = 5;
    structSimulationParameters.numDiffractionPoints = 2000;
    structSimulationParameters.approximateImpulseResponseCutDB = -90;
    structSimulationParameters.approximateImpulseResponse = 0;
    structSimulationParameters.ditherRaytracing = 1;
    structSimulationParameters.speedOfSound = 1500;

    %% Now calculate the whole setup

    rotationObjectAzimuth = -90:2:90;
    structSensor.position = [ 0.3 0 0];
    numTargetRotations = length( rotationObjectAzimuth );

    dataStorageMatrix = {};

    PB = ProgressBar( numTargetRotations, 'Running Processing', 'cli');
    for cntPoseObject = 1 : numTargetRotations
      
        structMeshPreparation.orientation = [ 0 0 rotationObjectAzimuth(cntPoseObject) ];
        structSurface = prepareMeshSurface( structMeshPreparation, 0 );
        structSimulationResult = calculateImpulseResponseFast( structSensor, structSurface, structSimulationParameters );

        structSimulatorOutput = struct();
        structSimulatorOutput.impulseResponse = structSimulationResult.impulseResponse;
        structSimulatorOutput.impulseResponseDiffraction = structSimulationResult.impulseResponseDiffraction;
        structSimulatorOutput.impulseResponseRaytracing = structSimulationResult.impulseResponseRaytracing;
        structSimulatorOutput.sensorInfo = structSensor;

        structSimulatorOutput.curObjectPose = [ 0 0 rotationObjectAzimuth(cntPoseObject) ];
        structSimulatorOutput.allObjectPoses = rotationObjectAzimuth;
        dataStorageMatrix{ cntPoseObject } = structSimulatorOutput;


        count( PB);
    end
