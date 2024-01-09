function [ subresultDiffraction ] = calculateImpulseRespondeDiffraction(structSensor, structSurface, structSimulationParameters)


    micsTransformed = ( calcBeamRotation( structSensor.orientation,  structSensor.coordsReceivers' ) + structSensor.position(:) );
    emitterTransformed = ( calcBeamRotation( structSensor.orientation,  structSensor.coordsEmitter' ) + structSensor.position(:) )';
    
    importanceVertex = structSurface.BRDF( :, 1 );

    importanceVertex = importanceVertex - min(importanceVertex );
    importanceVertex = importanceVertex / max(importanceVertex(:));
    importanceVertex = normLin( importanceVertex.^4  );

    [ importanceVertexSort, indexSort ] = sort( importanceVertex);
    
    importanceVertexCDF = normLin( cumsum(importanceVertexSort + 0.00001*randn(size(importanceVertexSort)) ) );

    numRandomSamples = structSimulationParameters.numDiffractionPoints;
    samplesRandomUniform = rand( numRandomSamples, 1 );

    samplesRandomImportance = round( interp1( importanceVertexCDF, 1:length( importanceVertexCDF ), samplesRandomUniform ) );

    vertexDiffraction = structSurface.surfaceVertices( indexSort( samplesRandomImportance ), : );
    vertexDiffractionNormals = structSurface.surfaceNormals( indexSort( samplesRandomImportance ), : );
    vectorEmitterToDiffractionNormed = ( vertexDiffraction - emitterTransformed ) ./ vecnorm( vertexDiffraction -  emitterTransformed, 2, 2 );
    materialsDiffraction = structSurface.surfaceMaterial( indexSort( samplesRandomImportance ), : );
    conicalAngleIntersect = acosd( diag( vertexDiffractionNormals*vectorEmitterToDiffractionNormed' ) );

    idxAngleValid = find( conicalAngleIntersect > 90 & conicalAngleIntersect < 270 );
    nMics = size( structSensor.coordsReceivers, 1 );
    nFrequencies = length( structSimulationParameters.vecFreqSim );
    pointsReflectedDiffraction = vertexDiffraction( idxAngleValid, : );

    strengthsReflectedDiffractionMaterial = repmat( structSurface.structMeshPreparation.materialStrengthDiffractive, [ size(pointsReflectedDiffraction,1) 1 nMics] );
    nPointsDiffracted = size( pointsReflectedDiffraction, 1 );
    pointsDistanceDiffraction = zeros( nPointsDiffracted, nMics );
    for cntMic = 1 : nMics
        pointsDistanceDiffraction( :, cntMic ) = sqrt( sum( (emitterTransformed - pointsReflectedDiffraction).^2, 2 ) ) + sqrt( sum( (micsTransformed(:, cntMic)' - pointsReflectedDiffraction).^2, 2 ) );
    end

    pathLossDiff = 1 ./ (pointsDistanceDiffraction.^2 );
    pathLossRepped = permute( repmat( pathLossDiff, [ 1 1 nFrequencies ]), [ 1 3 2 ]);

    rangeRepped = permute( repmat( pointsDistanceDiffraction, [ 1 1 nFrequencies ]), [ 1 3 2 ]);
    freqRepped = repmat( structSimulationParameters.vecFreqSim, [ size(pointsReflectedDiffraction,1) 1 nMics] );
    alphaAbsorption = 0.038 * (freqRepped / 1000) - 0.3;
    pathlossAbsorption = 10.^( -( alphaAbsorption .* rangeRepped ) / 20 );

    strengthsReflectedDiffraction = strengthsReflectedDiffractionMaterial .* pathLossRepped .* pathlossAbsorption;

    subresultDiffraction = struct();
    subresultDiffraction.pointsReflected = pointsReflectedDiffraction;
    subresultDiffraction.strengthsReflected = strengthsReflectedDiffraction;
    subresultDiffraction.distancesReflected = pointsDistanceDiffraction;
    

end

