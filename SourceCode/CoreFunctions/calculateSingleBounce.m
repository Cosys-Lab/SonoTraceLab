function structReturn = calculateSingleBounce( structSurface, structBounceParameters, doplot )
    % This function calculates a single ray intersection

    if ( ~isfield( structBounceParameters,'pEmitter') )
        structBounceParameters.pEmitter = structBounceParameters.intersectionValidPos';
    end

    if ( ~isfield( structBounceParameters,'directionRay') )
        structBounceParameters.directionRay = structBounceParameters.intersectionDirectionReflected;
    end
    
    % Vertex preparation
    vert1 = structSurface.surfaceVertices(structSurface.surfaceFaces(:,1),:);
    vert2 = structSurface.surfaceVertices(structSurface.surfaceFaces(:,2),:);
    vert3 = structSurface.surfaceVertices(structSurface.surfaceFaces(:,3),:);

    norms1 = structSurface.surfaceNormals(structSurface.surfaceFaces(:,1),:);
    norms2 = structSurface.surfaceNormals(structSurface.surfaceFaces(:,2),:);
    norms3 = structSurface.surfaceNormals(structSurface.surfaceFaces(:,3),:);

    % Intersection calculation
    [ intersectionGuesses, ~] = TriangleRayIntersection(structBounceParameters.pEmitter, structBounceParameters.directionRay', vert1, vert2, vert3, 'LineType','ray', 'border','inclusive','eps',1e-12, 'planeType','two sided' );

    idxVertIntersect = find( intersectionGuesses);

    pointsIntersect = [  vert1(idxVertIntersect,1), vert1(idxVertIntersect,2), vert1(idxVertIntersect,3); ...
                    vert2(idxVertIntersect,1), vert2(idxVertIntersect,2), vert2(idxVertIntersect,3); ...
                    vert3(idxVertIntersect,1), vert3(idxVertIntersect,2), vert3(idxVertIntersect,3); ];

    normsIntersect = [  norms1(idxVertIntersect,1), norms2(idxVertIntersect,2), norms3(idxVertIntersect,3) ;...
                    norms1(idxVertIntersect,1), norms2(idxVertIntersect,2), norms3(idxVertIntersect,3) ;...
                    norms1(idxVertIntersect,1), norms2(idxVertIntersect,2), norms3(idxVertIntersect,3); ];

    BRFDIntersect = [ structSurface.BRDF( idxVertIntersect ) ; structSurface.BRDF( idxVertIntersect );  structSurface.BRDF( idxVertIntersect ) ]; 
    
    % Calculate angle with normal
    conicalAngleIntersect = acosd(normsIntersect * structBounceParameters.directionRay);
    idxAngleValid = find( conicalAngleIntersect > 90 & conicalAngleIntersect < 270 );
    
    if( isempty( idxAngleValid ) )
        isLastReflection = 1;
        rayDidHit = 0;
    else
        isLastReflection = 0;
        rayDidHit = 1;
    end

    % HERE YOU NEED TO ADD THE STUFF TO MAKE THE MULTIPLE REFLECTION
    % INTERSECTION PROBLEM GO AWAY. SELECT THE POINT WITH THE MINIMAL
    % DISTANCE TO THE SOURCE, AFTER YOU SELECTED FOR THE ONES WITH THE
    % VALID CONICAL ANGLES.
    
%     if( length( idxVertIntersect )>0 )
%         print('lets go');
%         % Here we know we found valid points from the point of view of the
%         % conical angle. Now, select the point that is closest:
%         pointsIntersectValidAngle = 
%         distancesPointsToSource = sum(( pointsIntersect - structBounceParameters.pEmitter ).^2,2);
%         [ ~, idxPointFinal ] = min( distancesPointsToSource );
%         
%     end
% try
    if( ~isempty(pointsIntersect) )
        distancesPointsToSource = sum(( pointsIntersect - structBounceParameters.pEmitter(:)' ).^2,2);
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
% catch e 
%     print('error')
% end
    % Return struct generation
%     structReturn = struct();
%     structReturn.intersectionValidPos = toVertVec(mean( pointsIntersect( idxAngleValid, : ) ));
%     structReturn.intersectionValidAngle = toVertVec( mean( conicalAngleIntersect( idxAngleValid ) ) );
%     structReturn.intersectionValidRange = toVertVec( sqrt( sum( ( structBounceParameters.pEmitter(:) - structReturn.intersectionValidPos(:) ).^2 ) ) );
%     structReturn.intersectionValidNormal = toVertVec( mean( normsIntersect( idxAngleValid, : ) ) );
%     structReturn.intersectionDirectionReflected = toVertVec( -2*(structReturn.intersectionValidNormal*structBounceParameters.directionRay') * structReturn.intersectionValidNormal + structBounceParameters.directionRay );
%     structReturn.intersectionBRDF = toVertVec( mean( BRFDIntersect( idxAngleValid, : ) ) );
%     structReturn.isLastReflection = isLastReflection;
%     structReturn.bounceNumber =  structBounceParameters.bounceNumber + 1;
%     structReturn.rayDidHit =  rayDidHit;

    structReturn = struct();
    structReturn.intersectionValidPos = toVertVec(mean( pointsIntersect( idxSortValid, : ) ));
    structReturn.intersectionValidAngle = toVertVec( mean( conicalAngleIntersect( idxSortValid ) ) );
    structReturn.intersectionValidRange = toVertVec( sqrt( sum( ( structBounceParameters.pEmitter(:) - structReturn.intersectionValidPos(:) ).^2 ) ) );
    structReturn.intersectionValidNormal = toVertVec( mean( normsIntersect( idxSortValid, : ) ) );
    structReturn.intersectionDirectionReflected = toVertVec( -2*(structReturn.intersectionValidNormal*structBounceParameters.directionRay') * structReturn.intersectionValidNormal + structBounceParameters.directionRay );
    structReturn.intersectionBRDF = toVertVec( mean( BRFDIntersect( idxSortValid, : ) ) );
    structReturn.isLastReflection = isLastReflection;
    structReturn.bounceNumber =  structBounceParameters.bounceNumber + 1;
    structReturn.rayDidHit =  rayDidHit;



    if( doplot == 1 )
        intersectionValidPos = structReturn.intersectionValidPos;
        intersectionValidNormal = structReturn.intersectionValidNormal;
        intersectionDirectionReflected = structReturn.intersectionDirectionReflected;
        hold on; plot3( intersectionValidPos(1), intersectionValidPos(2), intersectionValidPos(3), 'black.', 'markersize', 10); hold off
        
        vecPlot = [intersectionValidPos intersectionValidPos + 0.1*intersectionValidNormal];        
        hold on; plot3( vecPlot(1,:), vecPlot(2,:), vecPlot(3, :), 'black' ); hold off;

        vecPlot = [ structBounceParameters.pEmitter(:) intersectionValidPos ];
        hold on; plot3( vecPlot(1,:), vecPlot(2,:), vecPlot(3, :), 'r' ); hold off;

        vecPlot = [intersectionValidPos intersectionValidPos + 0.1*intersectionDirectionReflected];        
        hold on; plot3( vecPlot(1,:), vecPlot(2,:), vecPlot(3, :), 'b' ); hold off;    
    end

end