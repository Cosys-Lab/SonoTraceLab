function [ structSimulationResult ] = calculateImpulseResponseFast( structSensor, structSurface, structSimulationParameters )

    % This function calculates the impulse response for a certain
    % raytracing setup. It is the main function to call
    % This one attempts to accelerate the process
    % with a GPU-implementation of Möller and Trumbore's triangle ray 
    % intersection algorithm that works well when many directions are
    % calculated at once instead of one by one. 
  
    % These two lines fix some undefinedness on the sizes of position and orientation
    structSensor.orientation = structSensor.orientation(:);
    structSensor.position = structSensor.position(:)';


    % Do the raytracing part:
    [ subresultRaytracing ] = calculateImpulseResponeRaytracing(structSensor, structSurface, structSimulationParameters);

    % Do the Diffraction Part:
    [ subresultDiffraction ] = calculateImpulseRespondeDiffraction(structSensor, structSurface, structSimulationParameters);

    impulseResponseRaytracing = synthetizeImpulseResponseFFT( structSensor, structSimulationParameters, subresultRaytracing );
    impulseResponseDiffraction = synthetizeImpulseResponseFFT( structSensor, structSimulationParameters, subresultDiffraction );
    impulseResponse = impulseResponseRaytracing + impulseResponseDiffraction;
 

    pointsReflected = [ subresultRaytracing.pointsReflected ; subresultDiffraction.pointsReflected ];
    strengthsReflected = [ subresultRaytracing.strengthsReflected ; subresultDiffraction.strengthsReflected  ];
    distancesReflected = [ subresultRaytracing.distancesReflected ; subresultDiffraction.distancesReflected ];


    structSimulationResult = struct();
    structSimulationResult.impulseResponse = impulseResponse;
    structSimulationResult.impulseResponseDiffraction = impulseResponseDiffraction;
    structSimulationResult.impulseResponseRaytracing = impulseResponseRaytracing;
    structSimulationResult.subresultRaytracing = subresultRaytracing;
    structSimulationResult.subresultDiffraction = subresultDiffraction;
    structSimulationResult.pointsReflected = pointsReflected;
    structSimulationResult.strengthsReflected = strengthsReflected;
    structSimulationResult.distancesReflected = distancesReflected;

end


