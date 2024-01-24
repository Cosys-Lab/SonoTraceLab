function [ structSimulationResult ] = calculateImpulseResponse( structSensor, structSurface, structSimulationParameters )

    % This function calculates the impulse response for a certain
    % raytracing setup. It is the main function to call
  

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
    
    % Plot the directions on unit sphere
    if( structSimulationParameters.doPlot )
        figure; plot3( pointsSphere(1,:), pointsSphere(2,:), pointsSphere(3,:), '.');
            axis equal
            hold on;
                plot3( pointsSphereSelected(1,:), pointsSphereSelected(2,:), pointsSphereSelected(3,:), '.');
            hold off
            axis equal
            grid on;
            xlabel( 'X-axis' );
            ylabel( 'Y-axis' );
            zlabel( 'Z-axis' );
    end

    % Transform to degrees
    [azVec, elVec, ~] = cart2sph(pointsSphereSelected(1, :), pointsSphereSelected(2, :), pointsSphereSelected(3, :));
    azVecAzEl = rad2deg(azVec);
    elVecAzEl = rad2deg(elVec);    

   % This is the main loop, which is the most expensive
    PB = ProgressBar( length( elVecAzEl ), 'Computing Raytracing', 'cli');
    dataMask = nan( length( elVecAzEl ), 2, structSensor.nMics );
    lastReflectionPoints = nan( length( elVecAzEl ), 3 );
    dataNumBounces = nan( length( elVecAzEl ), 1 );
    parfor cntDirection = 1 : length( elVecAzEl )
        
        % Choose a direction to emit a ray
        azRayEmit = azVecAzEl( cntDirection );
        elRayEmit = elVecAzEl( cntDirection );
        
%         azRayEmit = 0;
%         elRayEmit = 0;        
%         
        % Add the direction of emission to the orientation
        emissionBeamDirection =  structSensor.orientation + [ 0 elRayEmit azRayEmit ]';
        directionRay = calcBeamRotation( emissionBeamDirection, [ 1 0 0]' );
        
        % Setup the struct for the single ray trace
        structBounceParameters = struct();
        structBounceParameters.pEmitter = emitterTransformed;
        structBounceParameters.directionRay = [];
        structBounceParameters.pReceiver = micsTransformed;
        structBounceParameters.bounceNumber = 0;
        structBounceParameters.MAXBOUNCES = 5;
        structBounceParameters.directionRay = directionRay;

        % Calculate the raypath:
        rayReturnStruct = calculateRaypath( structSurface, structBounceParameters,  structSimulationParameters.doPlot  );

        % If we did hit something, store the last reflection location, the
        % distance travelled and the strength
        if( rayReturnStruct.didHit )
            dataMask( cntDirection, :, : ) =  [ rayReturnStruct.distanceTravelled ; rayReturnStruct.reflectionStrength ];
            lastReflectionPoints( cntDirection, : ) = rayReturnStruct.lastValidBounce.intersectionValidPos;
            dataNumBounces(cntDirection) = rayReturnStruct.numBounces;
        end

        count(PB)
    end
    
    % Linear impulse response generation
    
    impulseResponse = zeros( structSimulationParameters.numSamplesImpresp, structSensor.nMics );
    sampleRate =  structSimulationParameters.sampleRateImpresp;
    idxHits = find( ~isnan( dataMask( :,1, 1 ) ) );
    numHits = length( idxHits );
    for cnthits = 1 : numHits
        curDistance = squeeze( dataMask( idxHits(cnthits), 1, : ) );
        curTime = curDistance / structSimulationParameters.speedOfSound;
        curSample = round( curTime * sampleRate );
        for cntMic = 1 : structSensor.nMics
            impulseResponse( curSample( cntMic ), cntMic ) = impulseResponse( curSample( cntMic ), cntMic ) + dataMask( idxHits(cnthits), 2, cntMic );
        end
    end

    % Some variables for plotting the reflections that are strong enough. A
    % bit cumbersome, I know.
    pointsReflected = lastReflectionPoints( idxHits, : );
    strengthsReflected = squeeze( dataMask( idxHits, 2, 1 ) );
    idxStrongReflections = find(strengthsReflected>0.1 );
    pointsReflectedOK = pointsReflected( idxStrongReflections, : );

    structSimulationResult = struct();
    structSimulationResult.impulseResponse = impulseResponse;
    structSimulationResult.azVecAzEl = azVecAzEl;
    structSimulationResult.elVecAzEl = elVecAzEl;
    structSimulationResult.dataMask = dataMask;
    structSimulationResult.numBounces = dataNumBounces;
    structSimulationResult.reflectionPoints.idxHits = idxHits;
    structSimulationResult.reflectionPoints.pointsReflected = pointsReflected;
    structSimulationResult.reflectionPoints.strengthsReflected = strengthsReflected;


    if( structSimulationParameters.doPlot )
        figure()       
            patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 1 0 0], 'EdgeAlpha', 0.3); 
            axis equal;
            xlabel( 'X-Axis' )
            ylabel( 'Y-Axis' )
            zlabel( 'Z-Axis' )
            grid on
            axis equal
            hold on;
                drawTriad( structSensor.position(:), structSensor.orientation(:), 0.05)
                plot3(  micsTransformed(1,:), micsTransformed(2,:), micsTransformed(3,:),'.')
                plot3( pointsReflectedOK(:,1), pointsReflectedOK(:,2), pointsReflectedOK(:,3), 'b.')
            hold off;


        figure()       
            patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 1 0 0], 'EdgeAlpha', 0.3); 
            axis equal;
            xlabel( 'X-Axis' )
            ylabel( 'Y-Axis' )
            zlabel( 'Z-Axis' )
            grid on
            axis equal
            hold on;
                plot3( pointsReflectedOK(:,1), pointsReflectedOK(:,2), pointsReflectedOK(:,3), 'b.')
            hold off;

    end
    
end


