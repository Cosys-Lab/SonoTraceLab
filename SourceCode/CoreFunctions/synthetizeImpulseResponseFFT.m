function [ impulseResponse ] = synthetizeImpulseResponseFFT( structSensor, structSimulationParameters, subResult )

    % This is an approach for generating the impulse response for each microphone using a frequency domain approach. The idea is that we generate for each reflection
    % and for each microphone an impulse response using FFTs, based on the strength of the received echo for each frequency. We use a
    % gaussian window to smooth the frequency response a bit.

    pointsReflected = [ subResult.pointsReflected ];
    strengthsReflected = [ subResult.strengthsReflected ];
    distancesReflected = [ subResult.distancesReflected ];

    sampleRate = structSimulationParameters.sampleRateImpresp; 
    numSamplesIRFilter = structSimulationParameters.numSamplesIRFilter;
    freqVecIRFilter = (0:numSamplesIRFilter-1) * (sampleRate/numSamplesIRFilter);
    vecFreqSim = structSimulationParameters.vecFreqSim;
    IRGaussfiltAlpha = structSimulationParameters.IRFilterGaussAlpha;
    windowsGaussianFilter = gausswin( numSamplesIRFilter, IRGaussfiltAlpha );
    IRFilterFreqDomPhase = exp( -1i*linspace( 0, 0*pi, numSamplesIRFilter ));
           
    % Linear impulse response generation    
    impulseResponse = zeros( structSimulationParameters.numSamplesImpresp, structSensor.nMics );
    sampleRate =  structSimulationParameters.sampleRateImpresp;
    numHits = size( pointsReflected, 1 );
    
    % Approximation calculation
    dbCutApprox = structSimulationParameters.approximateImpulseResponseCutDB;
    energyReflectedNormed = normLog( sum( sum( strengthsReflected.^2, 2 ) , 3 ), dbCutApprox );
    idxsPointsApprox = find( energyReflectedNormed > (dbCutApprox+1) );
    numHitsApprox = size( idxsPointsApprox, 1 );

    if( structSimulationParameters.approximateImpulseResponse == 1 )
        for cntHitsApprox = 1 : numHitsApprox
            curDistance = distancesReflected( idxsPointsApprox(cntHitsApprox), : );
            curStrengths = squeeze( strengthsReflected( idxsPointsApprox(cntHitsApprox), :, : ) );
            curTime = curDistance / structSimulationParameters.speedOfSound;
            curSample = round( curTime * sampleRate );
    
            IRFilterFreqDomMagnitude = interp1( vecFreqSim, curStrengths, freqVecIRFilter, 'linear', 0 );
            IRFilterFreqDom = IRFilterFreqDomMagnitude .* IRFilterFreqDomPhase(:);
            IRFilterImpresp = fftshift( ifft( IRFilterFreqDom,[],1, 'symmetric'),1 );
            IRFilterImpresp = IRFilterImpresp .* windowsGaussianFilter;
            for cntMic = 1 : structSensor.nMics
                curSplStart = max( 0, curSample( cntMic ) - numSamplesIRFilter/2 );
                curSplStop = curSplStart + numSamplesIRFilter -1;
                impulseResponse( curSplStart:curSplStop , cntMic ) = impulseResponse( curSplStart:curSplStop , cntMic )  + IRFilterImpresp( :, cntMic );
            end
    
        end
    else
        for cntHits = 1 : numHits
            curDistance = distancesReflected( cntHits, : );
            curStrengths = squeeze( strengthsReflected( cntHits, :, : ) );
            curTime = curDistance / structSimulationParameters.speedOfSound;
            curSample = round( curTime * sampleRate );
    
            IRFilterFreqDomMagnitude = interp1( vecFreqSim, curStrengths, freqVecIRFilter, 'linear', 0 );
            IRFilterFreqDom = IRFilterFreqDomMagnitude .* IRFilterFreqDomPhase(:);
            IRFilterImpresp = fftshift( ifft( IRFilterFreqDom,[],1, 'symmetric'),1 );
            IRFilterImpresp = IRFilterImpresp .* windowsGaussianFilter;
            for cntMic = 1 : structSensor.nMics
                curSplStart = max( 0, curSample( cntMic ) - numSamplesIRFilter/2 );
                curSplStop = curSplStart + numSamplesIRFilter -1;
                impulseResponse( curSplStart:curSplStop , cntMic ) = impulseResponse( curSplStart:curSplStop , cntMic )  + IRFilterImpresp( :, cntMic );
            end
    
        end
    end
end

