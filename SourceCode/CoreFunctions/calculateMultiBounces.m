function structReturnBounces = calculateMultiBounces( structSurface, structBouncesParameters, doPlot )
    % This function calculates a set of rays

    % Vertex preparation
    nDirections = size(structBouncesParameters.directionRays, 2);
    vert1 = structSurface.surfaceVertices(structSurface.surfaceFaces(:,1),:);
    vert2 = structSurface.surfaceVertices(structSurface.surfaceFaces(:,2),:);
    vert3 = structSurface.surfaceVertices(structSurface.surfaceFaces(:,3),:);

    norms1 = structSurface.surfaceNormals(structSurface.surfaceFaces(:,1),:);
    norms2 = structSurface.surfaceNormals(structSurface.surfaceFaces(:,2),:);
    norms3 = structSurface.surfaceNormals(structSurface.surfaceFaces(:,3),:);

    % Intersection calculation
    nFaces = size(vert1, 1);
    indices = reshape(1 : (3 * nFaces), 3, nFaces)';
    interleavedVert = zeros(3*nFaces, 3);
    interleavedVert(indices(:), :) = [vert1; vert2; vert3];

    nRounds = floor(nDirections / structBouncesParameters.numberOfDirectionsPerCall);
    nDirectionsLast = mod(nDirections,structBouncesParameters.numberOfDirectionsPerCall);
    idxFaceIntersect_row = [];
    idxFaceIntersect_col = [];
    for cntRound = 0 : nRounds - 1
        curIndexes = cntRound * (structBouncesParameters.numberOfDirectionsPerCall) + 1 : (cntRound + 1) * (structBouncesParameters.numberOfDirectionsPerCall);
        curStartPoints = structBouncesParameters.startPoints(:, curIndexes );
        curDirectionRays = structBouncesParameters.directionRays(:, curIndexes );
        [curIdxVertIntersect_row, curIdxVertIntersect_col] = mollertrumborewrapper(interleavedVert', curStartPoints, curDirectionRays);
        curIdxVertIntersect_col = curIdxVertIntersect_col + (cntRound * structBouncesParameters.numberOfDirectionsPerCall);
        idxFaceIntersect_row = [idxFaceIntersect_row ; curIdxVertIntersect_row];
        idxFaceIntersect_col = [idxFaceIntersect_col ; curIdxVertIntersect_col];
    end
    lastIndexes = nDirections - nDirectionsLast + 1 :nDirections;
    lastStartPoints = structBouncesParameters.startPoints(:, lastIndexes );
    lastDirectionRays = structBouncesParameters.directionRays(:, lastIndexes );
    [lastIdxVertIntersect_row, lastIdxVertIntersect_col] = mollertrumborewrapper(interleavedVert', lastStartPoints, lastDirectionRays);
    lastIdxVertIntersect_col = lastIdxVertIntersect_col + (nDirections - nDirectionsLast);
    idxFaceIntersect_row = [idxFaceIntersect_row ; lastIdxVertIntersect_row];
    idxFaceIntersect_col = [idxFaceIntersect_col ; lastIdxVertIntersect_col];

    structReturnBounces = struct();
    structReturnBounces.bounceNumber = structBouncesParameters.bounceNumber + 1;
    % structReturnBounces.rayDidHit = zeros(nDirections, 1);
    % structReturnBounces.intersectionValidPos = zeros(3, nDirections);
    % structReturnBounces.intersectionDirectionReflected = zeros(3, nDirections);
    % structReturnBounces.intersectionValidNormal = zeros(3, nDirections);
    % allBounceStructCells = cell(nDirections, 1);

    tempRayDidHit = zeros(nDirections, 1);
    tempIntersectionValidPos = zeros(3, nDirections);
    tempIntersectionDirectionReflected = zeros(3, nDirections);
    tempIntersectionValidNormal = zeros(3, nDirections);
    allBounceStructCells = cell(nDirections, 1);

    parfor cntRay = 1 : nDirections
        idxFaceIntersect = idxFaceIntersect_row(idxFaceIntersect_col==cntRay);

        pointsIntersect = [ vert1(idxFaceIntersect,1), vert1(idxFaceIntersect,2), vert1(idxFaceIntersect,3); ...
                            vert2(idxFaceIntersect,1), vert2(idxFaceIntersect,2), vert2(idxFaceIntersect,3); ...
                            vert3(idxFaceIntersect,1), vert3(idxFaceIntersect,2), vert3(idxFaceIntersect,3); ];

        normsIntersect = [ norms1(idxFaceIntersect,1), norms2(idxFaceIntersect,2), norms3(idxFaceIntersect,3) ;...
                           norms1(idxFaceIntersect,1), norms2(idxFaceIntersect,2), norms3(idxFaceIntersect,3) ;...
                           norms1(idxFaceIntersect,1), norms2(idxFaceIntersect,2), norms3(idxFaceIntersect,3); ];

        % BRFDIntersect = [ structSurface.BRDF( idxVertIntersect, : ) ; structSurface.BRDF( idxVertIntersect, : );  structSurface.BRDF( idxVertIntersect, : ) ]; 
        % materialsIntersect = [ structSurface.surfaceMaterial( idxVertIntersect, : ) ; structSurface.surfaceMaterial( idxVertIntersect, : );  structSurface.surfaceMaterial( idxVertIntersect, : ) ];
        
        idxVertForMaterialAndBRDF = structSurface.surfaceFaces( idxFaceIntersect, 1 );
        BRFDIntersect = [ structSurface.BRDF(idxVertForMaterialAndBRDF, : ) ; structSurface.BRDF( idxVertForMaterialAndBRDF, : );  structSurface.BRDF( idxVertForMaterialAndBRDF, : ) ];
        materialsIntersect = [ structSurface.surfaceMaterial( idxVertForMaterialAndBRDF, : ) ; structSurface.surfaceMaterial( idxVertForMaterialAndBRDF, : );  structSurface.surfaceMaterial( idxVertForMaterialAndBRDF, : ) ];

        % Calculate angle with normal
        conicalAngleIntersect = acosd(normsIntersect * structBouncesParameters.directionRays(:, cntRay));
        idxAngleValid = find( conicalAngleIntersect > 90 & conicalAngleIntersect < 270 );

        if( isempty( idxAngleValid ) )
            rayDidHit = 0;
            curStructReturn = struct();
            curStructReturn.bounceNumber = structBouncesParameters.bounceNumber + 1;
            curStructReturn.rayDidHit = rayDidHit;
            allBounceStructCells{cntRay} = curStructReturn;
            tempRayDidHit(cntRay)= rayDidHit;
        else
            rayDidHit = 1;

            if( ~isempty(pointsIntersect) )
                distancesPointsToSource = sum(( pointsIntersect - structBouncesParameters.startPoints(:, cntRay)' ).^2,2);
                % This is a quick hack to make sure that the unvalid points (where
                % conical angle is not valid) have way larger distances in this
                % test function than the valid ones. However, this invalidates
                % distancesPointsToSource!!! they are not correct distances
                % anymore, and should not be used.
                distancesPointsToSource( idxAngleValid ) = distancesPointsToSource( idxAngleValid ) * 0.001;
                [ distSort, idxSort ] = sort( distancesPointsToSource );
                idxSortValid = idxSort(1:3);
            else
                idxSortValid = idxAngleValid;
            end
    
            curStructReturn = struct();
            curStructReturn.bounceNumber = structBouncesParameters.bounceNumber + 1;
            curStructReturn.intersectionValidPos = toVertVec(mean( pointsIntersect( idxSortValid, : ) ));
            curStructReturn.intersectionValidAngle = toVertVec( mean( conicalAngleIntersect( idxSortValid ) ) );
            curStructReturn.intersectionValidRange = toVertVec( sqrt( sum( ( structBouncesParameters.startPoints(:, cntRay) - curStructReturn.intersectionValidPos(:) ).^2 ) ) );
            curStructReturn.intersectionValidNormal = toVertVec( mean( normsIntersect( idxSortValid, : ) ) );
            curStructReturn.intersectionDirectionReflected = toVertVec( -2*(curStructReturn.intersectionValidNormal*structBouncesParameters.directionRays(:, cntRay)') * curStructReturn.intersectionValidNormal + structBouncesParameters.directionRays(:, cntRay) );
            
            % figure(1212);
            % hold on;
            % plot3(curStructReturn.intersectionValidPos(1), curStructReturn.intersectionValidPos(2), curStructReturn.intersectionValidPos(3), 'b.', 'markersize',100)
            % hold off
            % BRFDIntersect
            for cntFreq = 1 : size(BRFDIntersect, 2)
                curStructReturn.intersectionBRDF(:, cntFreq) = toVertVec( mean( BRFDIntersect( idxSortValid, cntFreq ) ) );
                curStructReturn.intersectionMaterial(:, cntFreq) = toVertVec( mean( materialsIntersect( idxSortValid, cntFreq ) ) );
            end     


            curStructReturn.rayDidHit = rayDidHit;
            allBounceStructCells{cntRay} = curStructReturn;
            tempRayDidHit(cntRay)= rayDidHit;
            tempIntersectionValidPos(:, cntRay) = curStructReturn.intersectionValidPos;
            tempIntersectionDirectionReflected(:, cntRay) = curStructReturn.intersectionDirectionReflected;
            tempIntersectionValidNormal(:, cntRay) = curStructReturn.intersectionValidNormal;

            
        end      
    end

    structReturnBounces.rayDidHit = tempRayDidHit;
    structReturnBounces.intersectionValidPos = tempIntersectionValidPos;
    structReturnBounces.intersectionDirectionReflected = tempIntersectionDirectionReflected;
    structReturnBounces.intersectionValidNormal = tempIntersectionValidNormal;
    structReturnBounces.allBounceStructCells = allBounceStructCells;

    if( doPlot == 1 )
        
        intersectionValidPos = structReturnBounces.intersectionValidPos;
        intersectionValidNormal = structReturnBounces.intersectionValidNormal;
        intersectionDirectionReflected = structReturnBounces.intersectionDirectionReflected;
        hold on; plot3( intersectionValidPos(1, :), intersectionValidPos(2, :), intersectionValidPos(3, :), 'green.', 'markersize', 30); hold off
        
        hold on; plot3([intersectionValidPos(1,:)' (intersectionValidPos(1,:)' + 0.1*intersectionValidNormal(1,:)')]',...
                       [intersectionValidPos(2,:)' (intersectionValidPos(2,:)' + 0.1*intersectionValidNormal(2,:)')]',...
                       [intersectionValidPos(3,:)' (intersectionValidPos(3,:)' + 0.1*intersectionValidNormal(3,:)')]', 'black' ); hold off;

        hold on; plot3([structBouncesParameters.startPoints(1,:)' intersectionValidPos(1,:)']',...
                       [structBouncesParameters.startPoints(2,:)' intersectionValidPos(2,:)']',...
                       [structBouncesParameters.startPoints(3,:)' intersectionValidPos(3,:)']', 'g' ); hold off;


        hold on; plot3([structBouncesParameters.startPoints(1,:)' (structBouncesParameters.startPoints(1,:)' + 0.1*structBouncesParameters.directionRays(1,:)')]',...
                       [structBouncesParameters.startPoints(2,:)' (structBouncesParameters.startPoints(2,:)' + 0.1*structBouncesParameters.directionRays(2,:)')]',...
                       [structBouncesParameters.startPoints(3,:)' (structBouncesParameters.startPoints(3,:)' + 0.1*structBouncesParameters.directionRays(3,:)')]', 'r' ); hold off;   

        hold on; plot3([intersectionValidPos(1,:)' (intersectionValidPos(1,:)' + 0.1*intersectionDirectionReflected(1,:)')]',...
                       [intersectionValidPos(2,:)' (intersectionValidPos(2,:)' + 0.1*intersectionDirectionReflected(2,:)')]',...
                       [intersectionValidPos(3,:)' (intersectionValidPos(3,:)' + 0.1*intersectionDirectionReflected(3,:)')]', 'b' ); hold off;    
    end

    
end