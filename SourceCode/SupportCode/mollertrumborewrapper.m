function [outIdxVertIntersect_row, outIdxVertIntersect_col]  = mollertrumborewrapper(interleavedVertices, startPoints, directionRays)    
    
lastIntersectionGuesses = mollertrumbore(single(interleavedVertices), single(startPoints), single(directionRays)); 

    [curIdxVertIntersect_row, curIdxVertIntersect_col] = find(lastIntersectionGuesses > 0);     
    inds = sub2ind( size( lastIntersectionGuesses ), curIdxVertIntersect_row, curIdxVertIntersect_col );
    intersectionsFoundCompressed = lastIntersectionGuesses(inds);    
    examplesFound = dec2bin( intersectionsFoundCompressed, 32 )- '0';    
    numRows = size( examplesFound, 1 );
    outIdxVertIntersect_row = [];
    outIdxVertIntersect_col = [];
    for cntRow = 1 : numRows
        curVertIdx = curIdxVertIntersect_row( cntRow );
        curExpansion = examplesFound( cntRow, : );
        idxRaysHit = find( curExpansion > 0);
        for cnt = 1 : length( idxRaysHit )            
            foundRayIdx = (curIdxVertIntersect_col(cntRow)-1)*32 + (33-idxRaysHit( cnt ));
            outIdxVertIntersect_row = [outIdxVertIntersect_row ; curVertIdx];
            outIdxVertIntersect_col = [outIdxVertIntersect_col ; foundRayIdx];
            % fprintf( 'Found one: IdxVert: %d, IdxRay: %d\n', curVertIdx-1, foundRayIdx-1 );            
        end        
    end

end