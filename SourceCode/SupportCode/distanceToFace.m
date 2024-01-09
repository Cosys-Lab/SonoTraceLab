function distance = distanceToFace(point, faceVertices)
    % Calculate the normal vector of the plane defined by the face vertices
    v1 = faceVertices(1, :);
    v2 = faceVertices(2, :);
    v3 = faceVertices(3, :);

    % Calculate the normal vector of the plane
    normalVector = cross(v2 - v1, v3 - v1);
    normalVector = normalVector / norm(normalVector);

    % Calculate the distance from the point to the plane
    distance = abs(dot(normalVector, point - v1));
end