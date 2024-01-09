function [idxVertIntersect_row, idxVertIntersect_col] = mollertrumborewrapper(interleavedVertices, startPoints, directionRays)    
    [idxVertIntersect_row, idxVertIntersect_col] = mollertrumbore(single(interleavedVertices), single(startPoints), single(directionRays)); 
    idxVertIntersect_row = single(idxVertIntersect_row);
    idxVertIntersect_col = single(idxVertIntersect_col);
end