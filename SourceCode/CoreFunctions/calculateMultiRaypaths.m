function raysReturnStruct = calculateMultiRaypaths( structSurface, structBouncesParameters, doPlot )
    % This funcitons calculates a full bounce of a ray - it is multibounce.

    % Calculate first bounces
    nInitialDirections = size(structBouncesParameters.directionRays, 2);
    structBouncesParameters.startPoints = repmat(structBouncesParameters.pEmitter', 1, nInitialDirections);
    structReturnBounces = calculateMultiBounces( structSurface, structBouncesParameters, doPlot );  

    rayPathtravelled = zeros(nInitialDirections, 1);
    lastValidBouncesStructCells = cell(nInitialDirections, 1);
    allValidBounces = cell(nInitialDirections, structBouncesParameters.MAXBOUNCES);
    previousFullList = 1:nInitialDirections;
    firstHitIndexes = find(structReturnBounces.rayDidHit);
    for cntBounce = 2 : structBouncesParameters.MAXBOUNCES
        previousIntersectionHits = find(structReturnBounces.rayDidHit);
        if ~isempty(previousIntersectionHits)
            indexesToFill = previousFullList(previousIntersectionHits);
            previousValidBouncesStruct = [structReturnBounces.allBounceStructCells{previousIntersectionHits}];
            lastValidBouncesStructCells(indexesToFill) = num2cell(previousValidBouncesStruct);
            allValidBounces(indexesToFill, structReturnBounces.bounceNumber) = num2cell(previousValidBouncesStruct);
            rayPathtravelled(indexesToFill) = rayPathtravelled(indexesToFill) + [previousValidBouncesStruct.intersectionValidRange]';
    
            newStructBouncesParameters = {};
            newStructBouncesParameters.startPoints = structReturnBounces.intersectionValidPos(:, previousIntersectionHits);
            newStructBouncesParameters.bounceNumber = structReturnBounces.bounceNumber;
            newStructBouncesParameters.numberOfDirectionsPerCall = structBouncesParameters.numberOfDirectionsPerCall;
            newStructBouncesParameters.directionRays = structReturnBounces.intersectionDirectionReflected(:, previousIntersectionHits);
            structReturnBounces = calculateMultiBounces( structSurface, newStructBouncesParameters, doPlot );  
            previousFullList = indexesToFill;
        end
    end

    numReceivers = size( structBouncesParameters.pReceiver, 2 );
    numFreq = size(structSurface.surfaceMaterial, 2);
    raysReturnStruct = struct();
    raysReturnStruct.didHit = zeros(nInitialDirections, 1);
    raysReturnStruct.didHit(firstHitIndexes) = 1;
    raysReturnStruct.reflectionStrength = nan(nInitialDirections, numFreq, numReceivers);
    raysReturnStruct.allValidBounces = allValidBounces;
    raysReturnStruct.lastValidBounces = lastValidBouncesStructCells;
    onlyHitsOfLastValidBouncesStructCells = [lastValidBouncesStructCells{firstHitIndexes}];
    raysReturnStruct.numBounces = nan(nInitialDirections, 1);
    if( ~isempty(firstHitIndexes) )
        raysReturnStruct.numBounces(firstHitIndexes) = [onlyHitsOfLastValidBouncesStructCells.bounceNumber]';
    end
        
    raysReturnStruct.distanceTravelled = nan(nInitialDirections, numReceivers);
    raysReturnStruct.lastReflectionPoints = nan(nInitialDirections, 3);
    
    for cntRay = 1 : nInitialDirections

        if raysReturnStruct.didHit(cntRay)
            
            % This is the array adaptation. The idea is that we calculate the
            % distance from the last bounce to each microphone. This is mostly
            % correct.
            lastValidBouncesStruct = lastValidBouncesStructCells{cntRay};
            distanceToReceiver = sqrt( sum( ( lastValidBouncesStruct.intersectionValidPos - structBouncesParameters.pReceiver ).^2 ) );
            raysReturnStruct.distanceTravelled(cntRay, :) = rayPathtravelled(cntRay) + distanceToReceiver;
            
            vecReceiverToReflection = ( lastValidBouncesStruct.intersectionValidPos - structBouncesParameters.pReceiver ) ./ vecnorm( lastValidBouncesStruct.intersectionValidPos - structBouncesParameters.pReceiver, 2 );
            raysReturnStruct.lastReflectionPoints(cntRay, :) = lastValidBouncesStruct.intersectionValidPos;

            % Calculate the BRDF function of the last bounce
            angleReflection = acosd( lastValidBouncesStruct.intersectionDirectionReflected(:)' * vecReceiverToReflection );
            curReflectionStrengthPathLoss = 1./(raysReturnStruct.distanceTravelled(cntRay, :).^2 );
            
            
            for cntFreq = 1 : numFreq
                curBRDFExponent = -1/(2*lastValidBouncesStruct.intersectionBRDF(cntFreq)^2);               
                curMaterialValue = lastValidBouncesStruct.intersectionMaterial(cntFreq);               

                curFreq =  structSurface.vecFreqSim( cntFreq );
                alphaAbsorption = 0.038 * ( curFreq / 1000 ) - 0.3;
                pathlossAbsorption = 10.^(-(alphaAbsorption*raysReturnStruct.distanceTravelled(cntRay, :) )/20);

                % Calculate the strength
                curReflectionStrengthBRDF = exp( curBRDFExponent .* (angleReflection - 180 ).^2 );             
                raysReturnStruct.reflectionStrength(cntRay, cntFreq, :) = curReflectionStrengthBRDF .* curReflectionStrengthPathLoss .* curMaterialValue .* pathlossAbsorption;              
            end            
        end   
    end
end