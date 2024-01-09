function valueMixed = mixSigmoidProperties( X, a, b, V1, V2 )

  
    valSigmoid1 = 1 ./ (1 + exp(-a * (X - b)));
    valSigmoid2 = 1 - valSigmoid1;

    valueMixed = valSigmoid1*V1 + valSigmoid2*V2;


    % X = 0:0.001:10;
    %     b = 2
    %     a = 8
    %     valSigmoid1 = 1 ./ (1 + exp(-a * (X - b)));
    %     figure; plot( X, valSigmoid1 )
    % 
    %     figure; hist( valSigmoid1, 100)
