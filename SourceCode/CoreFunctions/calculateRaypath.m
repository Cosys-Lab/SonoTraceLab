
function rayReturnStruct = calculateRaypath( structSurface, structBounceParameters, doPlot )
    % This funcitons calculates a full bounce of a ray - it is multibounce.

    % Calculate first bounce: 
    structReturnBounce = calculateSingleBounce( structSurface, structBounceParameters, doPlot);
    % bounceHistory(1) = structReturnBounce;
    bounceNumber = 0;
    rayPathtravelled = 0;

    rayReturnStruct = struct();
    
    % Check if we dit a hit:
    if( structReturnBounce.rayDidHit == 0 )
        % the code for not hitting
        rayReturnStruct.didHit = 0;
        rayReturnStruct.distanceTravelled = nan;
        rayReturnStruct.reflectionStrength = nan;
        rayReturnStruct.lastValidBounce = [];
    else
        % We hti something, so now, we continue bouncing the beam until we
        % found lastReflection
        lastValidBounceStruct = structReturnBounce;
        while( ( structReturnBounce.isLastReflection == 0 ) && ( bounceNumber < ( structBounceParameters.MAXBOUNCES -1) ) )
            if( structReturnBounce.rayDidHit == 1 )
                % Update the distance:
                rayPathtravelled = rayPathtravelled + structReturnBounce.intersectionValidRange;
                lastValidBounceStruct = structReturnBounce;
            end
            % Try a new bounce:
            structReturnBounce = calculateSingleBounce( structSurface, structReturnBounce, doPlot);
            bounceNumber = bounceNumber + 1;
        end
        
        % This is the array adaptation. The idea is that we calculate the
        % distance from the last bounce to each microphone. This is mostly
        % correct.
        numReceivers = size( structBounceParameters.pReceiver, 2 )    ;
        distanceToReceiver = sqrt( sum( ( lastValidBounceStruct.intersectionValidPos - structBounceParameters.pReceiver ).^2 ) );
        rayPathtravelled = rayPathtravelled + distanceToReceiver;
        
        vecReceiverToReflection = ( lastValidBounceStruct.intersectionValidPos - structBounceParameters.pReceiver ) ./ vecnorm( lastValidBounceStruct.intersectionValidPos - structBounceParameters.pReceiver, 2 );
        
        % Calculate the BRDF function of the last bounce
        curBRDFExponent = -1/(2*lastValidBounceStruct.intersectionBRDF^2);
        angleReflection = acosd( lastValidBounceStruct.intersectionDirectionReflected(:)' * vecReceiverToReflection );
  
        % Calculate the strength
        curReflectionStrengthBRDF = exp( curBRDFExponent .* (angleReflection - 180 ).^2 );
        curReflectionStrengthPathLoss = 1./(rayPathtravelled.^2 );
        curReflectionStrength = curReflectionStrengthBRDF .* curReflectionStrengthPathLoss;

        % Return struct generation
        rayReturnStruct.didHit = 1;
        rayReturnStruct.distanceTravelled = rayPathtravelled;
        rayReturnStruct.reflectionStrength = curReflectionStrength;
        rayReturnStruct.lastValidBounce = lastValidBounceStruct;
        rayReturnStruct.numBounces = min( structBounceParameters.MAXBOUNCES, bounceNumber );
    end