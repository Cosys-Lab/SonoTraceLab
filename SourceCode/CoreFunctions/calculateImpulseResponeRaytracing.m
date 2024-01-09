function [ subresultRaytracing ] = calculateImpulseResponeRaytracing(structSensor, structSurface, structSimulationParameters)

    % Transform the sensor to it's position and orientation
    micsTransformed = ( calcBeamRotation( structSensor.orientation,  structSensor.coordsReceivers' ) + structSensor.position(:) );
    emitterTransformed = ( calcBeamRotation( structSensor.orientation,  structSensor.coordsEmitter' ) + structSensor.position(:) )';


    % Calculate the directions within the limits that are chosen, and
    % distribute the number of points there.
    nEqPointsInPart = structSimulationParameters.numberOfDirections;
    
    % Calculate the number of points that need to be distrubuted onto the
    % sphere to have the desired number of points in the area of interest
    surfacePart = ( deg2rad( structSimulationParameters.limitsAzimuth(2) ) - deg2rad( structSimulationParameters.limitsAzimuth(1) ) ) * ( sind( structSimulationParameters.limitsElevation(2) ) - sind( structSimulationParameters.limitsElevation(1) ) );
    scaler = 4*pi / surfacePart; 
    nEqPointsSphere = round( nEqPointsInPart * scaler );

    % Calculate the points on a sphere, and select the range
    pointsSphere = eq_point_set(2, nEqPointsSphere);
    [ pointsAz, pointsEl, ~ ] = cart2sph( pointsSphere(1,:), pointsSphere(2,:), pointsSphere(3,:) );
    idxPointsOK = find( pointsAz > deg2rad( structSimulationParameters.limitsAzimuth(1) ) & ...
                        pointsAz < deg2rad( structSimulationParameters.limitsAzimuth(2) ) & ...
                        pointsEl > deg2rad( structSimulationParameters.limitsElevation(1) ) & ... 
                        pointsEl < deg2rad( structSimulationParameters.limitsElevation(2) ) );
    pointsSphereSelected = pointsSphere(:, idxPointsOK);

    
    % % Plot the directions on unit sphere
    % if( structSimulationParameters.doPlot )
    %         hold on;
    %         plot3( pointsSphere(1,:), pointsSphere(2,:), pointsSphere(3,:), '.');
    %         axis equal
    %         hold on;
    %             plot3( pointsSphereSelected(1,:), pointsSphereSelected(2,:), pointsSphereSelected(3,:), '.');
    %         hold off
    %         axis equal
    %         grid on;
    %         xlabel( 'X-axis' );
    %         ylabel( 'Y-axis' );
    %         zlabel( 'Z-axis' );
    % end

    % Transform to degrees
    [azVec, elVec, ~] = cart2sph(pointsSphereSelected(1, :), pointsSphereSelected(2, :), pointsSphereSelected(3, :));
    azVecAzEl = rad2deg(azVec);
    elVecAzEl = rad2deg(elVec);    

    if( structSimulationParameters.ditherRaytracing == 1 )
        solidAngleDirection = rad2deg( 4*pi / nEqPointsSphere ) ;
        maxPerturbationAngle =  2*pi * sqrt(solidAngleDirection / pi);

        azVecAzEl = azVecAzEl + maxPerturbationAngle * ( rand( size( azVecAzEl ) ) - 0.5 ) .* cosd( elVecAzEl );
        elVecAzEl = elVecAzEl + maxPerturbationAngle * ( rand( size( azVecAzEl ) ) - 0.5 );
    end
    
    % [ x, y, z ] = sph2cart( deg2rad( azVecAzEl ), deg2rad( elVecAzEl ), ones( size( azVecAzEl ) ) );
    % pointsSphereClean = [ x; y; z ];
    % [ x, y, z ] = sph2cart( deg2rad( azVecAzEl2 ), deg2rad( elVecAzEl2 ), ones( size( azVecAzEl2 ) ) );
    % pointsSphereDither = [ x; y; z ];    

    % % Plot the directions on unit sphere
    % figure;
    %     plot3( pointsSphereClean(1,:), pointsSphereClean(2,:), pointsSphereClean(3,:), '.');
    %     axis equal
    %     hold on;
    %         plot3( pointsSphereDither(1,:), pointsSphereDither(2,:), pointsSphereDither(3,:), '.');
    %     hold off
    %     axis equal
    %     grid on;
    %     xlabel( 'X-axis' );
    %     ylabel( 'Y-axis' );
    %     zlabel( 'Z-axis' );
    %     xlim([.9 1.1])
    % ylim([-0.1 0.1])
    % zlim([-0.1 0.1])




   % This is the main loop, which is the most expensive
    % PB = ProgressBar( length( elVecAzEl ), 'Computing Raytracing', 'cli');
    dataMask = nan( length( elVecAzEl ), 2, structSensor.nMics );
    lastReflectionPoints = nan( length( elVecAzEl ), 3 );
    dataNumBounces = nan( length( elVecAzEl ), 1 );

    % Loop all directions to emit a ray to initially save this into vectors
    % and start the datastruct for the current batch of rays to calculate
    nDirections = numel( elVecAzEl );
    directionRays = zeros(3, nDirections);
    for cntDirection = 1 : nDirections
        azRayEmit = azVecAzEl( cntDirection );
        elRayEmit = elVecAzEl( cntDirection );

        % azRayEmit = 10;
        % elRayEmit = 10;
        % Add the direction of emission to the orientation
        emissionBeamDirection =  structSensor.orientation + [ 0 elRayEmit azRayEmit ]';
        directionRays( :, cntDirection ) = calcBeamRotation( emissionBeamDirection, [ 1 0 0]' );
    end

    % Setup the struct for the multi ray trace
    structBouncesParameters = struct();
    structBouncesParameters.pEmitter = emitterTransformed;
    structBouncesParameters.pReceiver = micsTransformed;
    structBouncesParameters.MAXBOUNCES = 3;
    structBouncesParameters.bounceNumber = 0;
    structBouncesParameters.directionRays = directionRays;
    structBouncesParameters.numberOfDirectionsPerCall = structSimulationParameters.numberOfDirectionsPerCall;
    
    raysReturnStruct = calculateMultiRaypaths( structSurface, structBouncesParameters, structSimulationParameters.doPlot );
    % hitIndexes = find(1 : nDirections ==  raysReturnStruct.didHit);
    hitIndexes = find(1 ==  raysReturnStruct.didHit);
    dataMask(hitIndexes, 1, :) = raysReturnStruct.distanceTravelled(hitIndexes, :);

    lastReflectionPoints = raysReturnStruct.lastReflectionPoints;
    dataNumBounces = raysReturnStruct.numBounces;
   

    % Some variables for plotting the reflections that are strong enough. A
    % bit cumbersome, I know.
    pointsReflected = lastReflectionPoints( hitIndexes, : );
    strengthsReflected = squeeze( raysReturnStruct.reflectionStrength( hitIndexes, :, : ) );
    distancesReflected = squeeze( raysReturnStruct.distanceTravelled( hitIndexes, : ) );
    idxStrongReflections = find(strengthsReflected(:,1,1) > 0.1 );
    pointsReflectedOK = pointsReflected( idxStrongReflections, :, : );

    subresultRaytracing = struct();
    subresultRaytracing.pointsReflected = pointsReflected;
    subresultRaytracing.strengthsReflected = strengthsReflected;
    subresultRaytracing.distancesReflected = distancesReflected;
    subresultRaytracing.idxStrongReflections = idxStrongReflections;
    subresultRaytracing.pointsReflectedStrong = pointsReflectedOK;
    subresultRaytracing.azVecAzEl = azVecAzEl;
    subresultRaytracing.elVecAzEl = elVecAzEl;
    subresultRaytracing.dataMask = dataMask;
    subresultRaytracing.numBounces = dataNumBounces;
    subresultRaytracing.reflectionPoints.idxHits = hitIndexes;
    subresultRaytracing.reflectionPoints.pointsReflected = pointsReflected;
    subresultRaytracing.reflectionPoints.strengthsReflected = strengthsReflected;
end

